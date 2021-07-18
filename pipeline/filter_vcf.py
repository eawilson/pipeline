#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob
from collections import defaultdict
from bisect import bisect_right
try:
    from contextlib import nullcontext
except ImportError: # <= 3.6
    class nullcontext(object):
        def __init__(self, wrapped):
            self.wrapped = wrapped

        def __enter__(self):
            return self.wrapped

        def __exit__(self, exc_type, exc_value, exc_traceback):
            pass



def filter_vcf():
    """Cell free pipeline2 filter variants by bed file.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('input_vcf', help="Path of the input vcf file.")
    parser.add_argument("-b", "--bed", help="Path of bed file of regions to be included in the output.")
    parser.add_argument("-o", "--output", help="Path of output vcf.", default="-", required=False)
    args = parser.parse_args()

    if not args.bed.endswith(".bed"):
        args.bed = glob.glob(f"{args.bed}/*bed")
        if len(args.bed) != 1:
            sys.exit("{args.bed} does not contain an unambiguous bed file")
        args.bed = args.bed[0]

    startstops = defaultdict(list)
    with open(args.bed, "rt") as f_in:
        for row in f_in:
            row = row.strip().split("\t")
            if len(row) < 3:
                sys.exit("{args.bed} is not a valid bed file")
            startstops[row[0]].append((int(row[1]) + 1, int(row[2])))

    starts_by_contig = defaultdict(list)
    stops_by_contig = defaultdict(list)
    for contig in startstops.keys():
        merged = []
        for start, stop in sorted(startstops[contig]):
            if not merged or start > merged[-1][1] + 1:
                merged.append([start, stop])
            else:
                merged[-1][1] = max(merged[-1][1], stop)
        starts_by_contig[contig], stops_by_contig[contig] = zip(*merged)


    with (open(args.input_vcf, "rt") if args.input_vcf != "-" else nullcontext(sys.stdin)) as f_in:
        with (open(args.output, "wt") if args.output != "-" else nullcontext(sys.stdout)) as f_out:
            for row in f_in:
                if row.startswith("#"):
                    f_out.write(row)
                    continue

                fields = row.split("\t")
                contig = fields[0]
                starts = starts_by_contig.get(contig)
                if starts is None:
                    continue
                stops = stops_by_contig[contig]

                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                len_ref = len(ref)
                len_alt = len(alt)
                 
                if len_ref == 1 and len_alt > 1: # insertion
                    i = bisect_right(starts, pos)
                    if i > 0 and stops[i - 1] >= pos + 1:
                        f_out.write(row)

                else:
                    if len_ref > 1 and ref[0] == alt[0]: # deletion or deletion-insertion with padding
                        pos += 1
                        len_alt -= 1
                        len_ref -= 1

                    i = bisect_right(starts, pos + len_ref - 1)
                    if i > 0 and stops[1 - 1] >= pos:
                        f_out.write(row)
                        


def main():
    try:
        filter_vcf()
    except OSError as e:
        # File input/output error. This is not an unexpected error therfore
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

