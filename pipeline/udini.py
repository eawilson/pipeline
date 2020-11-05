import pdb
import argparse
import gzip
import json
import sys
from contextlib import closing
from itertools import chain
from collections import defaultdict

R1 = 0
R2 = 1

NAME = 0
SEQ = 1
PLUS = 2
QUAL = 3



def order_fastqs(fastqs):
    pairs = defaultdict(list)
    for fastq in fastqs:
        with fqopen(fastq, "rt") as f:
            qname = f.readline().split()[0]
            pairs[qname].append(fastq)
    return list(chain(*(pair for pair in pairs.values())))



def fqopen(fn, *args, **kwargs):
    return (gzip.open if fn.endswith(".gz") else open)(fn, *args, **kwargs)



def edit_distance(one, two):
    tot = 0
    for c1, c2 in zip(one, two):
        if c1 != c2:
            tot += 1
    return tot



def udini(input_fastqs,
          output_fastq="output.fastq",
          interleaved=False,
          umi=None,
          umi_length=0,
          umi_stem_length=0,
          umi_sequences="",
          statistics="stats.json",
          min_read_length=50,
          max_consecutive_ns=2,
          rtrim = 0):
    
    # reversed so we can pop the fastqs in the original order.
    input_fastqs = list(reversed(order_fastqs(input_fastqs)))
    
    if not output_fastq.endswith(".fastq") and \
       not output_fastq.endswith(".fastq.gz"):
        sys.exit("Output must be a .fastq or .fastq.gz file")
    
    if not(interleaved) and len(input_fastqs) % 2:
        sys.exit("Must be an even number of paired fastqs")

    if umi == "thruplex":
        umi_length = 6
        umi_stem_length = 11

    elif umi == "thruplex_hv":
        umi_length = 7
        umi_stem_length = 1
        umi_sequences = "AAGCTGA ACAACGA ACTCGTA ATTGCTC CGAGTAC CGCTAAT " \
                        "CTAGTAG GACATCG GTCTCTG TACCTCA TCTGGTA TGAACGG " \
                        "ACGACTC ATCTGGA CAATAGC CCTAGGT CGTCTCA GAGTCTC " \
                        "GGCAATG TCCACTA TCTCCAT TGTCAAC TGTGTCT TTGTAGT "
    
    elif umi == "prism":
        umi_length = 8
        umi_stem_length = 0
        umi_sequences = "GAGACGAT GCACAACT TTCCAAGG GCGTCATT CGCATGAT " \
                        "GAAGGAAG ACGGAACA ACTGAGGT CGGCTAAT TGAAGACG " \
                        "GCTATCCT GTTACGCA TGGACTCT AGCGTGTT ATCCAGAG " \
                        "GATCGAGT CTTAGGAC TTGCGAAG GTGCCATA CTGTTGAC " \
                        "TCGCTGTT GATGTGTG TTCGTTGG ACGTTCAG AAGCACTG " \
                        "TTGCAGAC GTCGAAGA CAATGTGG ACCACGAT ACGACTTG " \
                        "GATTACCG ACTAGGAG"
    elif umi is not None:
        sys.exit(f"'{umi}' is not a known UMI type")
    
    
    valid = set(umi for umi in umi_sequences.split())
    if any(len(umi) != umi_length for umi in valid):
        sys.exit(f"Not all UMI sequences are {umi_length} nt long")

    if min_read_length < 2 * (umi_length + umi_stem_length) + 1:
        min_read_length = 2 * (umi_length + umi_stem_length) + 1
    max_consecutive_ns = "N" * max_consecutive_ns
    
    total_reads = 0
    invalid_umi_reads = [0, 0]
    invalid_short_reads = 0
    invalid_n_reads = 0
    fastqs = [None, None]
    lines = [[None, None, None, None], [None, None, None, None]]
    with fqopen(output_fastq, "wt") as f_out:
        while input_fastqs:
            with closing(fqopen(input_fastqs.pop(), "rt")) as fastqs[R1]:
                with closing(fqopen(input_fastqs.pop(), "rt") if not interleaved else fastqs[R1]) as fastqs[R2]:
                    
                    eof = False
                    while not eof:
                        for i in range(8):
                            read = i // 4
                            line = fastqs[read].readline()
                            if not line:
                                eof = True
                                break
                            lines[read][i % 4] = line
                        
                        if eof:
                            break

                        total_reads += 1
                        if min_read_length and (len(lines[R1][SEQ]) < min_read_length or len(lines[R1][SEQ]) < min_read_length):
                            invalid_short_reads += 1
                            continue
                        
                        if max_consecutive_ns and (max_consecutive_ns in lines[read][SEQ] or max_consecutive_ns in lines[read][SEQ]):
                            invalid_n_reads += 1
                            continue

                        for read in (R1, R2):
                            lines[read][NAME] = lines[read][NAME].split(" ")[0].rstrip()

                        if umi_length:
                            invalid_umi = False
                            umis = ["", ""]
                            for read in (R1, R2):
                                umi = lines[read][SEQ][:umi_length]                            
                                if valid and umi not in valid:
                                    best = nextbest = umi_length
                                    for potential in valid:
                                        ed = edit_distance(potential, umi)
                                        if ed < best:
                                            nextbest = best
                                            best = ed
                                            corrected = potential
                                        elif ed < nextbest:
                                            nextbest = ed
                                    if best > 1 or nextbest < 3:
                                        invalid_umi_reads[read] += 1
                                        invalid_umi = True
                                    umi = corrected
                                umis[read] = umi
                            
                            if invalid_umi:
                                continue 
                                
                            tag = "RX:Z:{}-{}\tQX:Z:{} {}".format(umis[R1],
                                                                    umis[R2],
                                                                    lines[R1][QUAL][:umi_length],
                                                                    lines[R2][QUAL][:umi_length])
                            for read in (R1, R2):
                                lines[read][NAME] += f" {tag}\n"
                                length = len(lines[read][SEQ]) - rtrim
                                lines[read][SEQ] = lines[read][SEQ][umi_length + umi_stem_length:length]
                                lines[read][QUAL] = lines[read][QUAL][umi_length + umi_stem_length:length]
                                    
                        else:
                            for read in (R1, R2):
                                lines[read][NAME] = "{}\n".format(lines[read][NAME])
                                if rtrim:
                                    length = len(lines[read][SEQ]) - rtrim
                                    lines[read][SEQ] = lines[read][SEQ][:length]
                                    lines[read][QUAL] = lines[read][QUAL][:length]
                                    
                        if lines[R1][NAME] != lines[R2][NAME]:
                            sys.exit("Mismatched paired reads, names don't match")
                        f_out.writelines(chain(*lines))
                            
                    if i:
                        sys.exit("Truncated fastq")
    
    try:
        with open(statistics, "rt") as f:
            stats = json.load(f)
    except OSError:
        stats = {}
    stats["total_reads"] = total_reads
    stats["invalid_short_reads"] = invalid_short_reads
    stats["invalid_n_reads"] = invalid_n_reads
    if umi_sequences:
        stats["invalid_umi_reads_r1"] = invalid_umi_reads[R1]
        stats["invalid_umi_reads_r2"] = invalid_umi_reads[R2]
    with open(statistics, "wt") as f:
        json.dump(stats, f, sort_keys=True, indent=4)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Input fastq files.")
    parser.add_argument("-o", "--output", help="Output fastq file.", dest="output_fastq", default=argparse.SUPPRESS)
    parser.add_argument("-i", "--interleaved", help="Each input fastq contains alternating reads 1 and 2.", action="store_const", const=True, default=argparse.SUPPRESS)

    parser.add_argument("-u", "--umi", help="UMI type, allowed = thruplex, thruplex_hv, prism.", default=argparse.SUPPRESS)
    parser.add_argument("-l", "--umi-length", help="UMI length.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-k", "--umi-stem-length", help="UMI stem length.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-q", "--umi-sequences", help="UMI sequences.", default=argparse.SUPPRESS)
    parser.add_argument("-r", "--rtrim", help="Trim bases from the end of the read.", default=argparse.SUPPRESS)

    parser.add_argument("-m", "--min-read-length", help="Reads shoter than min-read-legth will be filtered.", default=argparse.SUPPRESS)
    parser.add_argument("-n", "--max-consecutive-ns", help="Reads containing more Ns than max-consecutive-ns will be filtered.", default=argparse.SUPPRESS)
    
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="statistics", default=argparse.SUPPRESS)
    
    args = parser.parse_args()
    try:
        udini(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

