import pdb
import argparse
import sys
import csv
from collections import defaultdict, Counter
from itertools import chain
from collections.abc import Sequence, Mapping
from math import sqrt

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
    parser.add_argument("-o", "--output", help="Output sam file, may be sam (default) or fastq.", default="-")
    parser.add_argument("-b", "--bed", help="Bed file of regions to include in output.")
    parser.add_argument("-m", "--max-fragment-size", help="Maximum fragment size to be considered a genuine read pair.", type=int, default=1000)
    args = parser.parse_args()
    
    if args.output.endswith(".sam") or args.output == "-":
        output_type = SAM
    elif args.output.endswith(".fastq"):
        output_type = FASTQ
    else:
        sys.exit(f"Invalid output file {args.output}, must be of type sam or fastq")
    
    unpaired = {}
    
    with (open(args.input_sam, "rt") if not args.output == "-" else nullcontext(sys.stdin)) as f_in:
        with (open(args.output, "wt") if not args.output == "-" else nullcontext(sys.stdout)) as f_out:
            for row in f_in:
                if row.startswith("@"):
                    if output_type == SAM:
                        f_out.write(row)
                    continue
                
                fields = row.split("\t")
                flag = int(fields[FLAG])
 
                # Secondary or supplementary read or both segments unmapped
                if flag & NON_PRIMARY or (flag & BOTH_UNMAPPED) == BOTH_UNMAPPED:
                    continue
                
                try:
                    mate = unpaired.pop(fields[QNAME])
                except KeyError:
                    unpaired[fields[QNAME]] = fields
                    continue
                



def main():
    try:
        filter_sam()
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

