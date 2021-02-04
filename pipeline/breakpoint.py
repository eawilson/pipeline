import pdb
import argparse
import json
import sys
import os
import glob
import csv
from collections import Counter
from itertools import chain
from pipeline.utils import string2cigar


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
READ1 = 0x40
READ2 = 0x80
SECONDARY = 0X100
FILTERED = 0x200
SUPPPLEMENTARY = 0x800
SEC_OR_SUP = SECONDARY | SUPPPLEMENTARY
BOTH_UNMAPPED = UNMAPPED | MATE_UNMAPPED

CONSUMES_REF = "MDN=X"
CONSUMES_READ = "MIS=X"

L = 0
R = 1



def bp(chrom, pos):
    return "{}:{}".format(chrom, pos)



def cigar_len(cig, ops):
    return sum(num for num, op in cig if op in ops)



def breakpoint(sam, output="translocations.tsv", mapq=10):

    breakpoints = Counter()
    with open(sam, "rt") as f_in:
        current_qname = ""
        read = []
        for row in chain(f_in, [""]):
            if row.startswith("@"):
                continue
            
            segment = row.split("\t")
            qname = segment[QNAME]
            if qname != current_qname:
                
                if read and len(read) < 2:
                    sys.exit("Sam file must be sorted by name.")
                
                read = [seg for seg in read if int(seg[MAPQ]) >= mapq]
                if len(set(seg[RNAME] for seg in read)) > 1:
                    for segment in read:
                        cigar = string2cigar(segment[CIGAR])
                        if cigar and (cigar[0][1] == "S") != (cigar[-1][1] == "S"):
                            pos = int(segment[POS])
                            if cigar[0][1] == "S":
                                breakpoints[bp(segment[RNAME], pos)] += 1
                            elif cigar[-1][1] == "S":
                                breakpoints[bp(segment[RNAME],  + pos + cigar_len(cigar, CONSUMES_REF))] += 1
                
                current_qname = qname
                read = []
            read.append(segment)
    
    with open(output, "wt") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["Location", "Reads"])
        for loc_count in sorted(breakpoints.items(), key=lambda x:x[1], reverse=True):
            writer.writerow(loc_count)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sam', help="Sam file, must be sorted by name.")
    parser.add_argument("-o", "--output", help="Output file.", default=argparse.SUPPRESS)
    parser.add_argument("-M", "--filter-mapq-less-than", help="Filter fragments with mapq less than.", dest="mapq", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        breakpoint(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

