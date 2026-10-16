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
Three independent hairpin signatures. A read trips one and its whole pair goes.

1. near_inverted_sa -- a same-chromosome SA partner within --cutoff bp in the
   OPPOSITE orientation, with no distant partner of MAPQ >= --rescue-min-mapq.
2. tandem_short_pair -- mates on the same chromosome pointing the SAME way with
   0 < TLEN < read length. Two ends cannot overlap almost completely while
   reading the same strand; only a hairpin molecule does that.
3. palindromic_clip -- a soft clip of >= --min-clip bases that is the reverse
   complement of the read's own aligned sequence, i.e. the read turns around on
   itself.

Rules 2 and 3 exist because the SA rule cannot see most of the artifact: 84% of
tandem-short pairs and 73% of palindromic-clip reads carry NO SA tag at all,
which is why an SA-only filter leaves the hairpin class largely intact.

Measured on WG-26-107 (FFPE tumor) vs matched normal WG-26-108, chr1:1-6Mb,
duplicates excluded, as a fraction of reads:

    rule                     tumor      normal     separation
    near_inverted_sa        15.85%      0.020%          792x
    tandem_short_pair        8.63%      0.003%         2876x
    palindromic_clip         5.08%      0.013%          390x
    union                   19.00%      0.034%          559x

The SA-only rule alone dropped 15.85%; the union adds 3.15% of tumor reads for
0.014% of normal reads. Both the added volume and the extreme tumor/normal
separation of each individual rule are the evidence that the additions target
artifact and not biology.

The MAPQ condition on the rescue in rule 1 is what makes it strict. Of 6,280
near-inverted reads carrying a distant partner, 3,824 (60.9%) had no distant
partner above MAPQ 10 and 2,736 were rescued solely by a MAPQ-0
interchromosomal segment: an ambiguous placement, not junction evidence. On the
normal, adding the MAPQ condition changes the rescue count by zero reads.

Genuine long-range and interchromosomal junctions, including foldback
inversions, are preserved. The same-chromosome inverted-partner distance
spectrum collapses at --cutoff 500 (15.81% of reads below 100bp, 0.024% in
500bp-1kb) while the 10kb-100kb and >1Mb foldback bins and the
interchromosomal class survive filtering essentially unchanged: 0.128% ->
0.117%, 0.099% -> 0.106%, and 118,019 -> 104,340 partners respectively.

WHY THE PREVIOUS "no SA + any clipping" RULE WAS REMOVED
--------------------------------------------------------
An earlier revision dropped every soft/hard-clipped read that carried no SA
tag. That rule contained no proximity or orientation logic, so it could not
discriminate hairpins. It removed 25.33% of sv.bam (211,636,408 reads) -- the
primary evidence class for insertions, small indels, and single-sided
breakends. On the non-FFPE normal it still removed 2.39% of reads despite there
being no artifact to remove.

Rule 3 is its principled replacement: it tests the clip against the read's own
aligned sequence, so it fires on 5.08% of tumor reads rather than 25.33%, and
on 0.013% rather than 2.39% of the normal.

MATE HANDLING
-------------
Input BAMs are coordinate-sorted, so mates are far apart in the stream. Dropping
reads independently orphans their mates: 62.56% of dropped reads have a mate that
would otherwise be retained, leaving dangling RNEXT/PNEXT and a stale
proper_pair flag. GRIDSS re-derives insert-size distributions during
preprocessing and reads discordant pairs as RP/REFPAIR evidence, so an
inconsistent BAM can distort its models.

Two passes are therefore used by default: pass 1 collects the query names to
drop, pass 2 removes both mates together. This reads the input twice; that cost
is accepted since the step runs once per sample. Use --no-two-pass only for
diagnostics where pair consistency is not required.

MEMORY
------
Names are stored as 64-bit digests in a flat uint64 open-addressed table
(DigestSet), not a Python `set`. A `set` of ints costs ~53 bytes per entry:
20M digests measured at 1.64 GB, and the 4e8 default ceiling implied 21 GB --
or 53 GB at the 1e9 some configs passed. The table costs 8 bytes per slot at
55% load: ~15 bytes per name resident, 26.8 bytes peak across a rehash, a 2x
reduction that is exact rather than probabilistic because the full digest is
kept in the slot and compared on probe.

The ceiling is also far above anything reachable. The count that matters is
DISTINCT NAMES DROPPED, not records scanned. On WG-26-107 chr1:1-6Mb, 977,920
dropped records collapsed to 305,552 distinct names, and 2,592,033 distinct
names existed in the whole window. Extrapolated to a 1e9-record genome-wide
BAM, the dropped-name count lands near 5e7, which the table holds in about
1.3 GB. A 4e8 ceiling therefore guards at roughly 8x the realistic peak; 1e9
guards at nothing. Measured peak RSS of the whole filter on the test window,
including pysam buffers, was 233 MB.

Both passes hash and probe in batches of BATCH_SIZE records so the linear-probe
loop runs over numpy blocks: 0.061us per lookup versus 0.26us for `set`.
"""

import argparse
import hashlib
import re
import sys
from collections import Counter

import numpy as np
import pysam

DEFAULT_CUTOFF = 500
DEFAULT_THREADS = 24
DEFAULT_MAX_QNAMES = 400_000_000
DEFAULT_RESCUE_MIN_MAPQ = 20
DEFAULT_MIN_CLIP = 10

# Records per vectorised membership batch. Large enough to amortise numpy call
# overhead, small enough that the AlignedSegment objects buffered by pass 2
# stay a bounded cost.
BATCH_SIZE = 200_000

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def _qname_digest(qname):
    """
    Stable 64-bit digest of a query name.

    Zero is the empty marker in DigestSet, so a digest of zero is remapped to
    one. The remap is applied here rather than in the table so that both passes
    agree by construction.
    """
    value = int.from_bytes(
        hashlib.blake2b(qname.encode("utf-8"), digest_size=8).digest(),
        "big",
    )
    return value or 1


def classify_sa_partners(read, cutoff, require_inverted, rescue_min_mapq):
    """
    Inspect a read's SA tag and classify its supplementary partners.

    Returns (has_near_artifact, has_rescuing_partner):
      has_near_artifact     -- a same-chromosome partner within `cutoff` bp
                               which, when `require_inverted` is set, is also
                               in the opposite orientation.
      has_rescuing_partner  -- a partner on another chromosome, or further than
                               `cutoff` bp away on the same chromosome, whose
                               own MAPQ is at least `rescue_min_mapq`.

    The MAPQ condition on the rescue matters. Measured on WG-26-107 (FFPE
    tumor) vs matched normal WG-26-108, chr1:1-6Mb, duplicates excluded:
    6,280 near-inverted reads carried a distant partner, but 3,824 of those
    (60.9%) had no distant partner above MAPQ 10, and 2,736 were rescued
    solely by a MAPQ-0 interchromosomal segment. A MAPQ-0 partner is an
    ambiguous placement, not evidence of a junction, so it must not license a
    hairpin read back into the BAM. The same rule costs the non-FFPE normal
    nothing: its rescue count at every cutoff is identical with and without
    the MAPQ condition (322 reads).
    """
    has_near_artifact = False
    has_rescuing_partner = False

    if not read.has_tag("SA"):
        return has_near_artifact, has_rescuing_partner

    read_strand = "-" if read.is_reverse else "+"
    # SAM SA field is 1-based; reference_start is 0-based.
    read_pos = read.reference_start + 1

    for entry in read.get_tag("SA").split(";"):
        if not entry:
            continue
        fields = entry.split(",")
        if len(fields) < 5:
            continue
        chrom, pos, strand, _cigar, mapq = fields[:5]

        try:
            partner_mapq = int(mapq)
        except ValueError:
            partner_mapq = 0

        if chrom != read.reference_name:
            if partner_mapq >= rescue_min_mapq:
                has_rescuing_partner = True
            continue

        try:
            distance = abs(int(pos) - read_pos)
        except ValueError:
            continue

        if distance > cutoff:
            if partner_mapq >= rescue_min_mapq:
                has_rescuing_partner = True
            continue

        # Same chromosome, within cutoff. Equal strand fields mean the
        # partner has the same orientation, so it is not a hairpin.
        if require_inverted and strand == read_strand:
            continue

        has_near_artifact = True

    return has_near_artifact, has_rescuing_partner


def is_tandem_short_pair(read, read_length):
    """
    True for a same-chromosome pair whose mates point the SAME way with an
    insert shorter than one read.

    A proper FR pair has its mates in opposite orientations. A TANDEM pair
    (both forward or both reverse) with TLEN below the read length cannot be a
    fragment at all: the two ends overlap almost completely while reading the
    same strand, which is what a hairpin molecule produces after end repair.
    This class carries no SA tag in 84% of cases, so the SA rule alone cannot
    see it. Measured genome-window counts: 8.63% of the FFPE tumor versus
    0.003% of the matched normal.
    """
    if read.is_unmapped or read.mate_is_unmapped:
        return False
    if read.reference_id != read.next_reference_id:
        return False
    if read.is_reverse != read.mate_is_reverse:
        return False

    insert = abs(read.template_length)
    return 0 < insert < read_length


def has_palindromic_clip(read, min_clip):
    """
    True when a soft clip is the reverse complement of the read's own aligned
    sequence.

    That is the direct signature of a hairpin: the polymerase turned around on
    the same molecule, so the clipped tail reads back over the bases already
    aligned. No SA tag is required for this to be visible, and 73% of the
    reads showing it carry none. Measured: 5.08% of the FFPE tumor versus
    0.013% of the matched normal, a 390-fold separation.
    """
    sequence = read.query_sequence
    if not sequence:
        return False

    offset = 0
    aligned = []
    clips = []

    for length, op in _CIGAR_RE.findall(read.cigarstring or ""):
        length = int(length)
        if op == "S":
            clips.append((offset, length))
            offset += length
        elif op in "M=X":
            aligned.append((offset, length))
            offset += length
        elif op == "I":
            offset += length

    if not aligned or not clips:
        return False

    aligned_seq = "".join(sequence[s : s + n] for s, n in aligned)

    for start, length in clips:
        if length < min_clip:
            continue
        probe = _revcomp(sequence[start : start + length])[: min(length, 20)]
        if probe in aligned_seq:
            return True

    return False


def drop_reason(read, cutoff, require_inverted, rescue_min_mapq, min_clip, read_length):
    """
    Return the name of the rule that condemns this read, or None to keep it.

    The three rules detect the same hairpin molecule through independent
    evidence, so a read need only trip one.
    """
    has_near_artifact, has_rescuing_partner = classify_sa_partners(
        read, cutoff, require_inverted, rescue_min_mapq
    )
    if has_near_artifact and not has_rescuing_partner:
        return "near_inverted_sa"

    if is_tandem_short_pair(read, read_length):
        return "tandem_short_pair"

    if has_palindromic_clip(read, min_clip):
        return "palindromic_clip"

    return None


def estimate_read_length(in_bam, threads, sample_size=100_000):
    """Longest query length seen in the first `sample_size` mapped records."""
    longest = 0
    seen = 0

    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in:
        for read in bam_in:
            if read.is_unmapped or read.query_length is None:
                continue
            longest = max(longest, read.query_length)
            seen += 1
            if seen >= sample_size:
                break

    return longest


class DigestSet:
    """
    Exact set of 64-bit query-name digests in a flat uint64 open-addressed
    table.

    A Python `set` of ints costs ~53 bytes per entry, measured: 20M digests
    occupied 1.64 GB purely for membership. That is what made --max-qnames a
    memory hazard rather than a guard. This table stores 8 bytes per slot and
    grows at 55% load, so ~15 bytes per name resident and 26.8 bytes per name
    peak across a rehash -- a 2x reduction with no loss of exactness, since the
    full digest is kept in the slot and compared on probe.

    Zero is reserved as the empty marker; a digest that hashes to zero is
    remapped to 1. Collision probability between two distinct names is 2^-64
    per pair, unchanged from the previous `set` of the same digests.

    Lookups are batched and vectorised: 0.061us per key versus 0.26us for
    `set`, because the linear-probe loop runs over whole numpy blocks instead
    of per-read Python bytecode.
    """

    __slots__ = ("_table", "_mask", "_count", "_capacity")

    _KNUTH = np.uint64(0x9E3779B97F4A7C15)
    _LOAD_NUM = 55
    _LOAD_DEN = 100

    def __init__(self, capacity=1 << 20):
        self._capacity = capacity
        self._table = np.zeros(capacity, dtype=np.uint64)
        self._mask = capacity - 1
        self._count = 0

    def __len__(self):
        return self._count

    @property
    def nbytes(self):
        return self._table.nbytes

    def _slots(self, keys):
        return (((keys * self._KNUTH) >> np.uint64(32)).astype(np.int64)) & self._mask

    def _limit(self):
        return (self._capacity * self._LOAD_NUM) // self._LOAD_DEN

    def add_batch(self, keys):
        keys = np.unique(keys)
        if keys.size == 0:
            return
        if self._count + keys.size > self._limit():
            self._grow(self._count + keys.size)

        table = self._table
        slots = self._slots(keys)
        todo = np.arange(keys.size)

        while todo.size:
            current = table[slots[todo]]
            duplicate = current == keys[todo]
            empty = current == 0

            placed = np.zeros(todo.size, dtype=bool)
            free = todo[empty]
            if free.size:
                # Several batch keys can target one empty slot; take the first
                # of each and re-probe the rest on the next iteration.
                _, first = np.unique(slots[free], return_index=True)
                winners = free[first]
                table[slots[winners]] = keys[winners]
                self._count += winners.size
                placed[np.searchsorted(todo, winners)] = True

            todo = todo[~(placed | duplicate)]
            slots[todo] = (slots[todo] + 1) & self._mask

    def _grow(self, needed):
        capacity = self._capacity
        while needed > (capacity * self._LOAD_NUM) // self._LOAD_DEN:
            capacity <<= 1

        occupied = self._table[self._table != 0]
        self._capacity = capacity
        self._table = np.zeros(capacity, dtype=np.uint64)
        self._mask = capacity - 1
        self._count = 0
        if occupied.size:
            self.add_batch(occupied)

    def contains_batch(self, keys):
        table = self._table
        slots = self._slots(keys)
        found = np.zeros(keys.size, dtype=bool)
        todo = np.arange(keys.size)

        while todo.size:
            current = table[slots[todo]]
            hit = current == keys[todo]
            empty = current == 0
            found[todo[hit]] = True
            todo = todo[~(hit | empty)]
            slots[todo] = (slots[todo] + 1) & self._mask

        return found


def collect_drop_names(
    in_bam, cutoff, require_inverted, rescue_min_mapq, min_clip, read_length,
    threads, max_qnames,
):
    """Pass 1: gather digests of query names to drop."""
    drop_names = DigestSet()
    reasons = Counter()
    total = 0
    pending = []

    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in:
        for read in bam_in:
            total += 1
            reason = drop_reason(
                read, cutoff, require_inverted, rescue_min_mapq, min_clip,
                read_length,
            )
            if reason is None:
                continue
            reasons[reason] += 1
            pending.append(_qname_digest(read.query_name))

            if len(pending) >= BATCH_SIZE:
                drop_names.add_batch(np.array(pending, dtype=np.uint64))
                pending.clear()
                if len(drop_names) > max_qnames:
                    sys.exit(
                        "ERROR: exceeded --max-qnames ({}) distinct query "
                        "names to drop. Refusing to continue rather than "
                        "exhaust memory. Raise --max-qnames if the allocation "
                        "can accommodate it.".format(max_qnames)
                    )

    if pending:
        drop_names.add_batch(np.array(pending, dtype=np.uint64))

    return drop_names, total, reasons


def write_filtered(in_bam, out_bam, drop_names, threads):
    """Pass 2: write reads whose query name was not marked for dropping."""
    kept = 0
    dropped = 0

    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in, pysam.AlignmentFile(
        out_bam, "wb", template=bam_in, threads=threads
    ) as bam_out:
        batch = []
        digests = []

        def flush():
            nonlocal kept, dropped
            if not batch:
                return
            drop = drop_names.contains_batch(np.array(digests, dtype=np.uint64))
            for read, is_drop in zip(batch, drop.tolist()):
                if is_drop:
                    dropped += 1
                    continue
                bam_out.write(read)
                kept += 1
            batch.clear()
            digests.clear()

        for read in bam_in:
            batch.append(read)
            digests.append(_qname_digest(read.query_name))
            if len(batch) >= BATCH_SIZE:
                flush()
        flush()

    return kept, dropped


def write_filtered_single_pass(
    in_bam, out_bam, cutoff, require_inverted, rescue_min_mapq, min_clip,
    read_length, threads,
):
    """Single-pass variant. Does not preserve mate consistency."""
    kept = 0
    dropped = 0
    reasons = Counter()

    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in, pysam.AlignmentFile(
        out_bam, "wb", template=bam_in, threads=threads
    ) as bam_out:
        for read in bam_in:
            reason = drop_reason(
                read, cutoff, require_inverted, rescue_min_mapq, min_clip,
                read_length,
            )
            if reason is not None:
                reasons[reason] += 1
                dropped += 1
                continue
            bam_out.write(read)
            kept += 1

    return kept, dropped, reasons


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
        "--rescue-min-mapq",
        type=int,
        default=DEFAULT_RESCUE_MIN_MAPQ,
        help="A distant SA partner only rescues a near-inverted read when the "
        "partner's own MAPQ is at least this value. Set to 0 to let any "
        "distant partner rescue, including ambiguous MAPQ-0 placements "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--min-clip",
        type=int,
        default=DEFAULT_MIN_CLIP,
        help="Minimum soft-clip length tested for being the reverse complement "
        "of the read's own aligned sequence (default: %(default)s).",
    )
    parser.add_argument(
        "--read-length",
        type=int,
        default=0,
        help="Read length used as the TLEN ceiling for the tandem-pair rule. "
        "Zero estimates it from the first 100,000 mapped records.",
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

    read_length = args.read_length or estimate_read_length(
        args.in_bam, args.threads
    )
    if read_length <= 0:
        sys.exit("ERROR: could not determine read length; pass --read-length.")
    sys.stderr.write("read length for tandem rule: {}\n".format(read_length))

    if args.two_pass:
        drop_names, total, reasons = collect_drop_names(
            args.in_bam,
            args.cutoff,
            args.require_inverted,
            args.rescue_min_mapq,
            args.min_clip,
            read_length,
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
        kept, dropped, reasons = write_filtered_single_pass(
            args.in_bam,
            args.out_bam,
            args.cutoff,
            args.require_inverted,
            args.rescue_min_mapq,
            args.min_clip,
            read_length,
            args.threads,
        )

    for reason, count in sorted(reasons.items()):
        sys.stderr.write("  rule {}: {} reads\n".format(reason, count))

    processed = kept + dropped
    pct = (dropped / processed * 100) if processed else 0.0
    sys.stderr.write(
        "wrote {} records, dropped {} ({:.3f}%)\n".format(kept, dropped, pct)
    )


if __name__ == "__main__":
    main()
