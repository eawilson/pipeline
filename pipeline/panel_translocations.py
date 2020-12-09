import pdb
import argparse
import json
import sys
import os
import glob
import csv
from collections import Counter, defaultdict
from itertools import chain
from covermi import bed, Gr, Entry
from pipeline import run


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



def t(rname1, rname2):
    if rname1 == rname2:
        return "-"
    rnames = [rname1[3:], rname2[3:]]
    try:
        return "t({};{})".format(*sorted(int(rname) for rname in rnames)) 
    except ValueError:
        return "t({};{})".format(*sorted(rnames)) 



def string2cigar(cigstr):
    if cigstr == "*":
        return ()
    
    cig = []
    num = ""
    for char in cigstr:
        if char.isnumeric():
            num += char
        else:
            try:
                cig.append((int(num), char))
            except ValueError:
                sys.exit(f"Malformed cigar string {cigstr}")
            num = ""
    if num:
        raise sys.exit(f"Malformed cigar string {cigstr}")
    return cig



def cigar_len(cig, ops):
    return sum(num for num, op in cig if op in ops)



def panel_translocations(sams, targets, statistics="stats.json", max_fragment_size=1000):

    if targets and os.path.isdir(targets):
        beds = glob.glob(f"{targets}/*.bed")
        if len(beds) == 0:
            sys.exit(f'"{targets}" does not contain a bedfile')
        elif len(beds) > 1:
            sys.exit(f'"{targets}" contains multiple bedfiles')
        targets = beds[0]
    
    targets = Gr(bed(targets))
    for target in targets:
        loc = f"{target.chrom}:{target.start}-{target.stop}"
        target.name = loc if target.name == "." else f"{target.name}_{loc}"
    
    
    with open("trans.tsv", "wt") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["Sample", "Bait", "Depth", "Translocation", "Primary VAF", "Supplementary VAF", "Primary Depth", "Supplementary Depth", "Ratio"])
        for sam in sams:
            print(sam)
            run(["samtools", "sort", "-n", "-o", "z.sam", sam])
            
            trans = defaultdict(lambda:defaultdict(Counter))
            depths = Counter()
            
            with open("z.sam", "rt") as f_in:
                current_qname = ""
                read = []
                for row in chain(f_in, [""]):
                    if row.startswith("@"):
                        continue
                    
                    segment = row.split("\t")
                    qname = segment[QNAME]
                    if qname != current_qname:
                        identify_translocations(read, targets, depths, trans, max_fragment_size)
                        current_qname = qname
                        read = []
                    read.append(segment)
            
            sample = sam.split(".")[0]
            written = False
            for bait, data in trans.items():
                for t, counts in data.items():
                    d = depths[bait]
                    r = min(counts["primary"], counts["supplementary"]) / max(counts["primary"], counts["supplementary"])
                    if r >= 0.2 and max(counts["primary"], counts["supplementary"]) >= 5:
                        writer.writerow([sample, bait, d, t, counts["primary"]/d, counts["supplementary"]/d, counts["primary"], counts["supplementary"], r])
                        written = True
            if not written:
                writer.writerow([sample])
    
    
    
    
    
    #try:
        #with open(statistics, "rt") as f:
            #old_stats = json.load(f)
    #except OSError:
        #old_stats = {}
    #old_stats.update(stats)
    #with open(statistics, "wt") as f:
        #json.dump(old_stats, f, sort_keys=True, indent=4)



def identify_translocations(read, targets, depths, trans, max_fragment_size):
    if not read:
        return

    primary = []
    supplementary = []
    for segment in read:
        segment[FLAG] = int(segment[FLAG])
        if not segment[FLAG] & SECONDARY: # Primary or supplementary read
            segment[CIGAR] = string2cigar(segment[CIGAR])
            segment[POS] = int(segment[POS])
            if segment[FLAG] & SUPPPLEMENTARY:
                supplementary.append(segment)
            else:
                primary.append(segment)
    
    if len(primary) != 2:
        sys.exit("Multiple primary reads in SAM file")
    
    if int(primary[L][MAPQ]) < 10 or int(primary[R][MAPQ]) < 10:
        return
    
    # Both mapped to same reference and pointing in opposite directions
    # therefore may be a concordant read pair.
    if not (primary[L][FLAG] & UNMAPPED) and not (primary[R][FLAG] & UNMAPPED) and \
        primary[L][RNAME] == primary[R][RNAME] and \
        primary[L][FLAG] & RC != primary[R][FLAG] & RC:
        
        if primary[L][FLAG] & RC:
            primary = primary[::-1]
        
        stop = primary[R][POS] + cigar_len(primary[R][CIGAR], CONSUMES_REF) - 1
        size = primary[R][POS] + cigar_len(primary[R][CIGAR], CONSUMES_READ) - primary[L][POS] - 1
        for num, op in primary[L][CIGAR]:
            if op in "IS":
                size += num
            elif op in "DN":
                size -= num
        if size < 0 or (max_fragment_size and size > max_fragment_size):
            size = 0
        
    else:
        size = 0
    
    if size:
        segments = [Entry(primary[L][RNAME], primary[L][POS], stop)]
        remaining = supplementary
    else:
        segments = []
        remaining = chain(primary, supplementary)
        
    for segment in remaining:
        start = int(segment[POS])
        segments.append(Entry(segment[RNAME], segment[POS], segment[POS] + cigar_len(segment[CIGAR], CONSUMES_REF) - 1))
    
    best_offset = None
    match = None
    for segment in segments:
        for target in targets.touched_by(segment):
            offset = abs(segment.start + segment.stop - target.start - target.stop)
            try:
                if offset > best_offset:
                    continue
            except TypeError:
                pass
            best_offset = offset
            match = target.name
            match_rname = segment.chrom
    
    if match:
        depths[match] += 1
        if primary[L][RNAME] != primary[R][RNAME]:
            trans[match][t(primary[L][RNAME], primary[R][RNAME])]["primary"] += 1
            
        for segment in supplementary:
            if segment[RNAME] != match_rname:
                trans[match][t(match_rname, segment[RNAME])]["supplementary"] += 1



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sams', nargs="+", help="Sam file, must be sorted by name.")
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="statistics", default=argparse.SUPPRESS)
    parser.add_argument("-b", "--bed", help="Bed file of on-target regions.", dest="targets", default=argparse.SUPPRESS, required=True)
    parser.add_argument("-m", "--max-fragment-size", help="Maximum fragment size to be considered a genuine read pair.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        panel_translocations(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

