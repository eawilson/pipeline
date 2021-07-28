import pdb
import argparse
import sys
import csv
import glob
from collections import defaultdict, Counter
from itertools import chain
from collections.abc import Sequence, Mapping
from math import sqrt
from bisect import bisect_right

from .utils import save_stats, string2cigar

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


QNAME = 0
FLAG = 1
RNAME = 2
POS = 3
MAPQ = 4
CIGAR = 5
RNEXT = 6
PNEXT = 7
TLEN = 8
SEQ = 9
QUAL = 10

UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
RC = 0x10
MATERC = 0x20
READ1 = 0x40
READ2 = 0x80
SECONDARY = 0X100
FILTERED = 0x200
SUPPPLEMENTARY = 0x800
NON_PRIMARY = SECONDARY | SUPPPLEMENTARY
BOTH_UNMAPPED = UNMAPPED | MATE_UNMAPPED
READXRCXUNMAPPEDX = READ1 | READ2 | RC | MATERC | UNMAPPED | MATE_UNMAPPED
READX = READ1 | READ2

CONSUMES_REF = "MDN=X"
CONSUMES_READ = "MIS=X"

L = 0
R = 1

RCOMPLEMENT = str.maketrans("ATGC", "TACG")
SAM = 0
FASTQ = 1


def cigar_len(cig, ops):
    return sum(num for num, op in cig if op in ops)    



def filter_sam():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', help="Input sam file.")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-b", "--bed", help="Bed file of regions to include in output.")
    args = parser.parse_args()
    
    if not args.bed.endswith(".bed"):
        args.bed = glob.glob(f"{args.bed}/*.bed")
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

    for contig, positions in list(startstops.items()):
        merged = []
        for start, stop in sorted(positions):
            if not merged or start > merged[-1][1] + 1:
                merged.append([start, stop])
            else:
                merged[-1][1] = max(merged[-1][1], stop)
        startstops[contig] = list(zip(*merged))
    
    
    chrom = None
    qnames = set()
    with (open(args.input_sam, "rt") if not args.input_sam == "-" else nullcontext(sys.stdin)) as f_in:
        for row in f_in:
            if row.startswith("@"):
                continue
            
            fields = row.split("\t")
            flag = int(fields[FLAG])
            
            if fields[RNAME] != chrom:
                chrom = fields[RNAME]
                #if chrom =="chr2":
                    #break
                print("Reading", chrom, file=sys.stderr)
            
            # Secondary or supplementary read or both segments unmapped
            if flag & NON_PRIMARY or (flag & BOTH_UNMAPPED) == BOTH_UNMAPPED:
                continue
            
            try:
                starts, stops = startstops.get(fields[RNAME])
            except TypeError:
                continue
            
            start = int(fields[POS])
            cigar = string2cigar(fields[CIGAR])
            stop = start + cigar_len(cigar, CONSUMES_REF) - 1
            
            i = bisect_right(starts, stop)
            if i > 0 and stops[i - 1] >= start:
                qnames.add(fields[QNAME])


    with open(args.output, "wt") as f_out:
        for qname in qnames:
            f_out.write(f"{qname}\n")



def main():
    try:
        filter_sam()
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

