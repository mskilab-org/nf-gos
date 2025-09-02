#!/usr/bin/env python3

import pysam
import sys
import re

in_bam = sys.argv[1]
out_bam = sys.argv[2]
bp_cutoff = 500
clip_re = re.compile(r"[0-9]+[SH]")  # match any soft/hard clip

def sa_close_to_primary(read, cutoff):
    """Return True if SA tag has same-chromosome coord within cutoff bp."""
    if not read.has_tag("SA"):
        return False
    for entry in read.get_tag("SA").split(';'):
        if not entry:
            continue
        chrom, pos, strand, cigar, mapq, nm = entry.split(',')
        if chrom == read.reference_name:
            if abs(int(pos) - (read.reference_start + 1)) <= cutoff:
                return True
    return False

with pysam.AlignmentFile(in_bam, "rb", threads=24) as bam_in, \
     pysam.AlignmentFile(out_bam, "wb", template=bam_in, threads=24) as bam_out:

    for read in bam_in:
        drop = False

        if read.has_tag("SA"):
            # Rule 2: drop if any SA entry is too close
            if sa_close_to_primary(read, bp_cutoff):
                drop = True
        else:
            # Rule 1: drop if no SA and has clipping
            if read.cigarstring and clip_re.search(read.cigarstring):
                drop = True

        if not drop:
            bam_out.write(read)