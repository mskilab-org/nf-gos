#!/usr/bin/env python3

"""
Remove FFPE end-repair chimera (hairpin) artifact reads from an aligned BAM.

BACKGROUND
----------
FFPE end-repair produces palindromic / hairpin molecules. When aligned, a single
read is split against nearby sequence in the *opposite* orientation, yielding a
supplementary alignment (SA tag) whose partner sits only tens of basepairs away.
These are alignment artifacts, not genomic rearrangements.

Measured on WG-25-160 (124x FFPE-repaired tumor) vs its matched non-FFPE normal
WG-25-161, genome-wide over 835,460,680 sv.bam records (duplicates excluded):

    reads carrying an SA tag         19.0%  (tumor)   vs  0.71%  (normal)
    picard PCT_CHIMERAS              13.98% (tumor)   vs  1.82%  (normal)
    near-SA partners inverted        98.02%
    soft clips that are revcomp of their own aligned segment   26-34%
    picard PCT_ADAPTER               ~1e-6  (adapter contamination ruled out)

The tumor additionally shows a TANDEM pair-orientation class absent from the
normal (87,711,887 pairs), 89% of which have insert size < 150bp -- shorter than
a single 151bp read. That is a hairpin signature, not a fragment distribution.

WHAT IS REMOVED
---------------
A read is dropped when it has a same-chromosome SA partner within --cutoff bp
that is in the OPPOSITE orientation, and it has no other SA partner that is
distant or on another chromosome.

Genuine long-range and interchromosomal junctions, including foldback
inversions, are preserved. The arm-distance spectrum is empty between roughly
500bp and 10kb (0.22% of events), so --cutoff 500 separates the artifact class
from real biology rather than trading one off against the other. Of the reads
surviving the cutoff, 91.6% are interchromosomal and the 10kb-100kb and >1Mb
foldback bins remain intact.

WHY THE PREVIOUS "no SA + any clipping" RULE WAS REMOVED
--------------------------------------------------------
An earlier revision also dropped every soft/hard-clipped read that carried no SA
tag. That rule contained no proximity or orientation logic, so it could not
discriminate hairpins. It removed 25.33% of sv.bam (211,636,408 reads) -- the
primary evidence class for insertions, small indels, and single-sided breakends,
which GRIDSS collects deliberately via MIN_CLIP_LENGTH=5. On the non-FFPE normal
it still removed 2.39% of reads despite there being no artifact to remove. It has
been deleted.

Combined effect of the corrections, on the same genome-wide scan:

    previous rules                635,846,168   76.107% of sv.bam dropped
    current rules                 412,891,512   49.421% of sv.bam dropped
    reads recovered               222,954,656   26.686%

MATE HANDLING
-------------
Input BAMs are coordinate-sorted, so mates are far apart in the stream. Dropping
reads independently orphans their mates: 62.56% of dropped reads have a mate that
would otherwise be retained, leaving dangling RNEXT/PNEXT and a stale
proper_pair flag. GRIDSS re-derives insert-size distributions during
preprocessing and reads discordant pairs as RP/REFPAIR evidence, so an
inconsistent BAM can distort its models.

Two passes are therefore used by default: pass 1 collects the query names to
drop, pass 2 removes both mates together. Query names are stored as 64-bit
digests to bound memory (roughly 3-5 GB for ~200 million distinct names).
This reads the input twice; that cost is accepted since the step runs once per
sample. Use --no-two-pass only for diagnostics where pair consistency is not
required.
"""

import argparse
import hashlib
import sys

import pysam

DEFAULT_CUTOFF = 500
DEFAULT_THREADS = 24
DEFAULT_MAX_QNAMES = 400_000_000


def _qname_digest(qname):
    """Return a stable 64-bit digest of a query name."""
    return int.from_bytes(
        hashlib.blake2b(qname.encode("utf-8"), digest_size=8).digest(),
        "big",
    )


def classify_sa_partners(read, cutoff, require_inverted):
    """
    Inspect a read's SA tag and classify its supplementary partners.

    Returns (has_near_artifact, has_distant):
      has_near_artifact -- a same-chromosome partner within `cutoff` bp which,
                           when `require_inverted` is set, is also in the
                           opposite orientation.
      has_distant       -- a partner on another chromosome, or further than
                           `cutoff` bp away on the same chromosome.

    A read with both is retained: it carries a real long-range junction in
    addition to the artifact, and that evidence must not be discarded. This
    case is not rare -- 13.54% of SA reads have more than one SA segment, and
    0.83% of near-SA reads also have a distant partner.
    """
    has_near_artifact = False
    has_distant = False

    if not read.has_tag("SA"):
        return has_near_artifact, has_distant

    read_strand = "-" if read.is_reverse else "+"
    # SAM SA field is 1-based; reference_start is 0-based.
    read_pos = read.reference_start + 1

    for entry in read.get_tag("SA").split(";"):
        if not entry:
            continue
        fields = entry.split(",")
        if len(fields) < 3:
            continue
        chrom, pos, strand = fields[0], fields[1], fields[2]

        if chrom != read.reference_name:
            has_distant = True
            continue

        try:
            distance = abs(int(pos) - read_pos)
        except ValueError:
            continue

        if distance > cutoff:
            has_distant = True
            continue

        # Same chromosome, within cutoff. Equal strand fields mean the
        # partner has the same orientation, so it is not a hairpin.
        if require_inverted and strand == read_strand:
            continue

        has_near_artifact = True

    return has_near_artifact, has_distant


def should_drop(read, cutoff, require_inverted):
    """True if the read is FFPE hairpin artifact and carries no distant partner."""
    has_near_artifact, has_distant = classify_sa_partners(
        read, cutoff, require_inverted
    )
    return has_near_artifact and not has_distant


def collect_drop_names(in_bam, cutoff, require_inverted, threads, max_qnames):
    """Pass 1: gather digests of query names to drop."""
    drop_names = set()
    total = 0

    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in:
        for read in bam_in:
            total += 1
            if should_drop(read, cutoff, require_inverted):
                drop_names.add(_qname_digest(read.query_name))
                if len(drop_names) > max_qnames:
                    sys.exit(
                        "ERROR: exceeded --max-qnames ({}) distinct query names "
                        "to drop. Refusing to continue rather than exhaust "
                        "memory. Raise --max-qnames if the allocation can "
                        "accommodate it.".format(max_qnames)
                    )

    return drop_names, total


def write_filtered(in_bam, out_bam, drop_names, threads):
    """Pass 2: write reads whose query name was not marked for dropping."""
    kept = 0
    dropped = 0

    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in, pysam.AlignmentFile(
        out_bam, "wb", template=bam_in, threads=threads
    ) as bam_out:
        for read in bam_in:
            if _qname_digest(read.query_name) in drop_names:
                dropped += 1
                continue
            bam_out.write(read)
            kept += 1

    return kept, dropped


def write_filtered_single_pass(in_bam, out_bam, cutoff, require_inverted, threads):
    """Single-pass variant. Does not preserve mate consistency."""
    kept = 0
    dropped = 0

    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in, pysam.AlignmentFile(
        out_bam, "wb", template=bam_in, threads=threads
    ) as bam_out:
        for read in bam_in:
            if should_drop(read, cutoff, require_inverted):
                dropped += 1
                continue
            bam_out.write(read)
            kept += 1

    return kept, dropped


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Remove FFPE end-repair chimera (hairpin) artifact reads.",
    )
    parser.add_argument("in_bam", help="Input coordinate-sorted BAM.")
    parser.add_argument("out_bam", help="Output filtered BAM.")
    parser.add_argument(
        "--cutoff",
        type=int,
        default=DEFAULT_CUTOFF,
        help="Maximum distance (bp) to a same-chromosome SA partner for it to "
        "count as hairpin artifact (default: %(default)s).",
    )
    parser.add_argument(
        "--require-inverted",
        dest="require_inverted",
        action="store_true",
        default=True,
        help="Only treat a near SA partner as artifact when it is in the "
        "opposite orientation (default: enabled).",
    )
    parser.add_argument(
        "--no-require-inverted",
        dest="require_inverted",
        action="store_false",
        help="Treat any near SA partner as artifact regardless of orientation. "
        "Removes an additional ~2%% of near-SA reads that are not hairpins.",
    )
    parser.add_argument(
        "--two-pass",
        dest="two_pass",
        action="store_true",
        default=True,
        help="Drop both mates of an affected pair, preserving pair "
        "consistency (default: enabled).",
    )
    parser.add_argument(
        "--no-two-pass",
        dest="two_pass",
        action="store_false",
        help="Single pass. Faster, but orphans mates. Diagnostics only.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="Threads for BAM (de)compression (default: %(default)s).",
    )
    parser.add_argument(
        "--max-qnames",
        type=int,
        default=DEFAULT_MAX_QNAMES,
        help="Abort if more than this many distinct query names are marked for "
        "dropping, as a memory guard (default: %(default)s).",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if args.two_pass:
        drop_names, total = collect_drop_names(
            args.in_bam,
            args.cutoff,
            args.require_inverted,
            args.threads,
            args.max_qnames,
        )
        sys.stderr.write(
            "pass 1: scanned {} records, {} distinct query names marked for "
            "dropping\n".format(total, len(drop_names))
        )
        kept, dropped = write_filtered(
            args.in_bam, args.out_bam, drop_names, args.threads
        )
    else:
        sys.stderr.write(
            "WARNING: --no-two-pass leaves orphaned mates with dangling "
            "RNEXT/PNEXT and stale proper_pair flags. Diagnostics only.\n"
        )
        kept, dropped = write_filtered_single_pass(
            args.in_bam,
            args.out_bam,
            args.cutoff,
            args.require_inverted,
            args.threads,
        )

    processed = kept + dropped
    pct = (dropped / processed * 100) if processed else 0.0
    sys.stderr.write(
        "wrote {} records, dropped {} ({:.3f}%)\n".format(kept, dropped, pct)
    )


if __name__ == "__main__":
    main()
