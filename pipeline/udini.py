import pdb
import argparse
import gzip

from contextlib import closing
from itertools import chain
import json

R1 = 0
R2 = 1

NAME = 0
SEQ = 1
QUAL = 3



def fqopen(fn, *args, **kwargs):
    return (gzip.open if fn.endswith(".gz") else open)(fn, *args, **kwargs)



def edit_distance(one, two):
    tot = 0
    for c1, c2 in zip(one, two):
        if c1 != c2:
            tot += 1
    return tot



def udini(input_fastqs, output_fastq="output.fastq", interleaved=False, replace=False, umi=None, umi_length=0, umi_stem_length=0, umi_sequences=""):
    input_fastqs = reversed(input_fastqs)
    
    if not output_fastq.endswith() and not output_fastq.endswith():
        raise ValueError("Output must be a .fastq or .fastq.gz file.")
    
    if umi == "thruplex":
        umi_length = 6
        umi_stem_length = 11

    elif umi == "thruplex_hv":
        umi_length = 7
        umi_stem_length = 1
        valid = set(["AAGCTGA",
                     "ACAACGA",
                     "ACTCGTA",
                     "ATTGCTC",
                     "CGAGTAC",
                     "CGCTAAT",
                     "CTAGTAG",
                     "GACATCG",
                     "GTCTCTG",
                     "TACCTCA",
                     "TCTGGTA",
                     "TGAACGG",
                     "ACGACTC",
                     "ATCTGGA",
                     "CAATAGC",
                     "CCTAGGT",
                     "CGTCTCA",
                     "GAGTCTC",
                     "GGCAATG",
                     "TCCACTA",
                     "TCTCCAT",
                     "TGTCAAC",
                     "TGTGTCT",
                     "TTGTAGT"])
    
    elif umi == "prism":
        umi_length = 8
        umi_stem_length = 0
        valid = set(["GAGACGAT", "GCACAACT",
                     "TTCCAAGG", "GCGTCATT",
                     "CGCATGAT", "GAAGGAAG",
                     "ACGGAACA", "ACTGAGGT",
                     "CGGCTAAT", "TGAAGACG",
                     "GCTATCCT", "GTTACGCA",
                     "TGGACTCT", "AGCGTGTT",
                     "ATCCAGAG", "GATCGAGT",
                     "CTTAGGAC", "TTGCGAAG",
                     "GTGCCATA", "CTGTTGAC",
                     "TCGCTGTT", "GATGTGTG",
                     "TTCGTTGG", "ACGTTCAG",
                     "AAGCACTG", "TTGCAGAC",
                     "GTCGAAGA", "CAATGTGG",
                     "ACCACGAT", "ACGACTTG",
                     "GATTACCG", "ACTAGGAG"])
        
    elif umi is None:
        valid = set(umi for umi in umi_sequences.split())
        if any(len(umi) != umi_length for umi in valid):
            raise ValueError(f"Not all UMI sequences are {umi_length} nt long")
    
    else:
        raise ValueError(f"'{umi}' is not a known UMI type")
    
    total_reads = 0
    invalid_reads = 0
    fastqs = [None, None]
    lines = [[None, None, None, None], [None, None, None, None]]
    with fqopen(output_fastq, "wt") as f_out:
        while input_fastqs:
            with closing(fqopen(input_fastqs.pop())) as fastqs[R1]:
                if interleaved and not input_fastqs:
                    raise ValueError("Must be an even number of paired fastqs")
                with closing(fqopen(input_fastqs.pop()) if not interleaved else fastqs[R1]) as fastqs[R2]:
                    
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
                        if umi_length:
                            invalid_umi = False
                            umis = ["", ""]
                            for read in (R1, R2):
                                umi = lines[read][SEQ][:umi_length]                            
                                if valid and umi not in valid:
                                    best = nextbest = len(umi)
                                    for potential in valid:
                                        ed = edit_distance(potential, umi)
                                        if ed <= best:
                                            nextbest = best
                                            best = ed
                                            corrected = potential
                                    if best > 1 or nextbest < 3:
                                        invalid_umi = True
                                        invalid_reads += 1
                                        break
                                    umi = corrected
                                umis[read] = umi
                            
                            if invalid_umi:
                               continue 
                                
                            if replace:
                                for read in (R1, R2):
                                    lines[read][SEQ] = umis[read] + lines[read][SEQ][umi_length:]
                            
                            else:
                                tag = "RX:Z:{}-{}\tQX:Z:{} {}".format(umis[R1],
                                                                      umis[R2],
                                                                      line[R1][QUAL][:umi_length],
                                                                      line[R2][QUAL][:umi_length])
                                for read in (R1, R2):
                                    lines[read][NAME] = "{} {}\n".format(lines[read][NAME].split(" ")[0], tag)
                                    lines[read][SEQ] = lines[read][SEQ][umi_length + umi_stem_length:]
                                    lines[read][QUAL] = lines[read][QUAL][umi_length + umi_stem_length:]
                                    
                        f_out.writelines(chain(*lines))
                            
                    if i:
                        raise RuntimeError("Truncated fastq")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Input fastq files.")
    parser.add_argument("-o", "--output", help="Output fastq file.", dest="output_fastq", default=argparse.SUPPRESS)
    parser.add_argument("-i", "--interleaved", help="Each input fastq contains alternating reads 1 and 2.", action="store_const", const=True, default=argparse.SUPPRESS)
    parser.add_argument("-r", "--replace", help="Replace UMI back into read after correction.", action="store_const", const=True, default=argparse.SUPPRESS)

    parser.add_argument("-u", "--umi", help="UMI type, allowed = thruplex, thruplex_hv, prism.", type=int, default=argparse.SUPPRESS)
    
    parser.add_argument("-l", "--umi-length", help="UMI length.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-k", "--umi-stem-length", help="UMI stem length.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-q", "--umi-sequences", help="UMI sequences.", default=argparse.SUPPRESS)
    
    parser.add_argument("-s", "--stats", help="Statistics file.", default=argparse.SUPPRESS) # NOT CODED YET
    
    args = parser.parse_args()
    udini(**vars(args))



if __name__ == "__main__":
    main()

