#!/usr/bin/env python3

"""
Annotate GRIDSS breakend records with FFPE hairpin-chimera geometry and apply an
advisory FILTER tag.

WHY THIS EXISTS SEPARATELY FROM THE SUPPORT-BASED FILTER
--------------------------------------------------------
Depth and mapping quality do not separate FFPE hairpin artifacts from real
structural variants. Measured on WG-25-160 (124x FFPE-repaired tumor):

    recurrent artifact junctions with >=3 independent fragments   87.7%
    split reads supporting artifacts with MAPQ >= 20              93.7%

So the artifacts are reproducible, confidently mapped and multiply supported.
Any threshold on SR/RP/QUAL will pass them while cutting genuine low-support
SVs. The matched normal provides no help via germline subtraction either, since
it was not FFPE-repaired and carries essentially no artifact background
(73 opposite-strand local split reads vs 277,732 in the tumor over the same
2Mb window).

What does separate them is junction geometry. The artifact is a hairpin
turnaround, so it is bounded in span and inverted in orientation:

    opposite-strand split reads with partner < 200bp    93.7%
    arm distance 500bp - 10kb                            0.22%   (empty valley)
    near-SA partners that are inverted                  98.02%   (genome-wide)

Real foldback inversions place their arms kilobases to megabases apart. A
threshold in the empty valley therefore separates the two classes instead of
trading sensitivity for specificity.

GRIDSS encodes junction geometry only in the ALT breakend string, for example
"CTTCTTTC[1:54783[". bcftools filter expressions cannot parse that, which is
why this is a script rather than an expression.

WHAT IS ADDED
-------------
INFO/FFPE_ARMDIST    Distance in bp between the two breakend arms, or -1 when
                     the partner is on another chromosome.
INFO/FFPE_INVERTED   Flag, set when the ALT bracket orientation indicates an
                     inverted (foldback-like) junction.
FILTER/FFPE_GEOM_CHIMERA
                     Set when the junction is same-chromosome AND arm distance
                     is below --max-armdist AND inverted AND unsupported in the
                     matched normal.

This is advisory. No record is removed. Interchromosomal junctions and
long-range foldbacks are never tagged.

CAVEAT
------
The 500bp boundary is derived from alignment topology on a single FFPE sample.
It has not yet been validated against copy-number breakpoints, which is the
strongest available orthogonal check (a real foldback should coincide with a
copy-number step; a hairpin artifact should not). Treat the tag as a triage aid
until that validation is done.
"""

import argparse
import re
import sys

import pysam

DEFAULT_MAX_ARMDIST = 500

# GRIDSS/VCF breakend ALT forms:
#   t[chr:pos[   t]chr:pos]   [chr:pos[t   ]chr:pos]t
BND_RE = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")


def parse_breakend(alt):
    """
    Parse a VCF breakend ALT string.

    Returns (mate_chrom, mate_pos, bracket) or (None, None, None) when the ALT
    is a single-sided breakend (GRIDSS writes these with a '.' and they carry no
    mate coordinate).
    """
    match = BND_RE.search(alt)
    if not match:
        return None, None, None
    bracket = "[" if "[" in alt else "]"
    return match.group(1), int(match.group(2)), bracket


def is_inverted(alt, bracket):
    """
    Decide whether a breakend pair is inverted.

    In VCF breakend notation the bracket direction encodes the orientation of
    the joined segment. A non-inverted (deletion-like or duplication-like)
    junction places the mate sequence so that the bracket faces away from the
    reference base; an inverted junction has both breakends pointing the same
    way. Concretely, for ALT of the form t[p[ or ]p]t the junction preserves
    orientation, while t]p] and [p[t indicate inversion.
    """
    if bracket is None:
        return False
    leading_bracket = alt[0] in "[]"
    if bracket == "[":
        # [p[t is inverted; t[p[ is not.
        return leading_bracket
    # ]p]t is not inverted; t]p] is.
    return not leading_bracket


def normal_supports(record, normal_index):
    """
    True when the matched normal shows any split-read or read-pair support.

    Returns False when there is no normal sample, so that tumor-only runs are
    still eligible for tagging on geometry alone.
    """
    if normal_index is None:
        return False

    sample = record.samples[normal_index]
    for key in ("SR", "RP", "ASSR", "ASRP"):
        value = sample.get(key)
        if value is None:
            continue
        if isinstance(value, tuple):
            if any(v for v in value if v):
                return True
        elif value:
            return True
    return False


def add_headers(header):
    header.info.add(
        "FFPE_ARMDIST",
        1,
        "Integer",
        "Distance in bp between breakend arms; -1 if interchromosomal",
    )
    header.info.add(
        "FFPE_INVERTED",
        0,
        "Flag",
        "Breakend orientation indicates an inverted (foldback-like) junction",
    )
    header.filters.add(
        "FFPE_GEOM_CHIMERA",
        None,
        None,
        "Short-range inverted junction consistent with an FFPE end-repair "
        "hairpin artifact: same chromosome, arm distance below threshold, "
        "inverted orientation, and no support in the matched normal. Advisory "
        "only; long-range foldbacks and interchromosomal junctions are not "
        "tagged.",
    )


def resolve_normal_index(header, normal_id):
    if normal_id is None:
        return None
    samples = list(header.samples)
    if normal_id not in samples:
        sys.stderr.write(
            "WARNING: normal sample '{}' not present in VCF (samples: {}). "
            "Proceeding without matched-normal evidence.\n".format(
                normal_id, ",".join(samples) or "none"
            )
        )
        return None
    return samples.index(normal_id)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Add advisory FFPE hairpin-chimera geometry annotations to "
        "a GRIDSS breakend VCF.",
    )
    parser.add_argument("in_vcf", help="Input VCF/BCF (may be bgzipped).")
    parser.add_argument("out_vcf", help="Output bgzipped VCF.")
    parser.add_argument(
        "--max-armdist",
        type=int,
        default=DEFAULT_MAX_ARMDIST,
        help="Tag same-chromosome inverted junctions with arm distance below "
        "this many bp (default: %(default)s).",
    )
    parser.add_argument(
        "--normal-id",
        default=None,
        help="Sample name of the matched normal. When given, junctions with "
        "normal support are not tagged.",
    )
    args = parser.parse_args(argv)

    with pysam.VariantFile(args.in_vcf) as vcf_in:
        header = vcf_in.header.copy()
        add_headers(header)
        normal_index = resolve_normal_index(header, args.normal_id)

        total = 0
        tagged = 0
        inverted_count = 0

        with pysam.VariantFile(args.out_vcf, "wz", header=header) as vcf_out:
            for record in vcf_in:
                total += 1
                # translate() rewrites the record against the new header in
                # place and returns None, so keep using `record`.
                record.translate(header)
                new_record = record

                alt = new_record.alts[0] if new_record.alts else ""
                mate_chrom, mate_pos, bracket = parse_breakend(alt)

                if mate_chrom is None:
                    # Single-sided breakend: no mate coordinate, no geometry.
                    vcf_out.write(new_record)
                    continue

                if mate_chrom != new_record.chrom:
                    armdist = -1
                else:
                    armdist = abs(mate_pos - new_record.pos)

                inverted = is_inverted(alt, bracket)

                new_record.info["FFPE_ARMDIST"] = armdist
                if inverted:
                    new_record.info["FFPE_INVERTED"] = True
                    inverted_count += 1

                if (
                    armdist >= 0
                    and armdist < args.max_armdist
                    and inverted
                    and not normal_supports(new_record, normal_index)
                ):
                    new_record.filter.add("FFPE_GEOM_CHIMERA")
                    tagged += 1

                vcf_out.write(new_record)

    sys.stderr.write(
        "records: {}  inverted: {}  tagged FFPE_GEOM_CHIMERA: {} ({:.3f}%)\n".format(
            total,
            inverted_count,
            tagged,
            (tagged / total * 100) if total else 0.0,
        )
    )


if __name__ == "__main__":
    main()
