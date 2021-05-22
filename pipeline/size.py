import pdb
import argparse
import sys
import os
from collections import Counter, defaultdict
from itertools import chain

from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages

from .utils import run, save_stats, string2cigar, CONSUMES_REF, CONSUMES_READ


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

LEFT = 0
RIGHT = 1



def cigar_len(cig, ops):
    return sum(num for num, op in cig if op in ops)



def do_sizing(input_sam,
              output_file="output.filtered.sam",
              stats_file="stats.json",
              max_fragment_size=1000,
              rnames="",
              name="",
              output=""):

    contigs = set(rnames.split())
    stats = {"total": Counter()}
    for rname in contigs:
        stats[rname] = Counter()

    with open(input_sam, "rt") as f_in:
        current_qname = ""
        read = []
        for row in f_in:
            if row.startswith("@"):
                continue
            
            segment = row.split("\t")
            qname = segment[QNAME]                    
            if qname != current_qname:
                size_read(read, stats, contigs, max_fragment_size)
                current_qname = qname
                read = []
            read.append(segment)
            
        size_read(read, stats, contigs, max_fragment_size)


    if not name:
        name = os.path.basename(input_sam).split(".")[0]
    if not output:
        output = f"{name}.sizes.pdf"
    
    median_fragment_size = {}
    with PdfPages(output) as pdf:
        for contig in chain(["total"], sorted(key for key in stats.keys() if key != "total")):
            if not stats[contig]:
                continue
            
            title = name if contig == "total" else f"{name} - {contig}"
            
            sizes, counts = zip(*sorted(stats[contig].items()))
            median_count = sum(counts) // 2
            for median, count in zip(sizes, counts):
                median_count -= count
                if median_count <= 0:
                    break
            median_fragment_size[contig] = median
        
            figure = Figure(figsize=(11.69,8.27))
            FigureCanvasPdf(figure)
            ax = figure.gca()
            
            ax.set_xlim(left=0, right=max_fragment_size or None)
            ax.plot(sizes, counts, "-", color="dodgerblue", linewidth=1)
            ax.get_yaxis().set_visible(False)
            ax.axvline(median, color="black", linewidth=0.5, linestyle=":")
            ax.set_xlabel("Fragment size (bp)", fontsize=10)
            ax.set_title(title, fontsize=12)
            
            pdf.savefig(figure)

    save_stats(stats_file, {"fragment_sizes": stats,
                            "median_fragment_size": median_fragment_size})



def size_read(read, stats, contigs, max_fragment_size):
    if not read:
        return ""
    
    primary = []
    for segment in read:
        segment[FLAG] = int(segment[FLAG])
        if not segment[FLAG] & SEC_OR_SUP: # Primary read
            primary.append(segment)

    if len(primary) != 2:
        sys.exit("SAM file not filtered by name or corrupt")
    
    size = 0
    # Both mapped to same reference and pointing in opposite directions
    # therefore may be a concordant read pair.
    if not (primary[0][FLAG] & UNMAPPED) and not (primary[1][FLAG] & UNMAPPED) and \
        primary[0][RNAME] == primary[1][RNAME] and \
        primary[0][FLAG] & RC != primary[1][FLAG] & RC:
        
        if primary[0][FLAG] & RC:
            primary = primary[::-1]
        
        right_cigar = string2cigar(primary[1][CIGAR])
        size = int(primary[1][POS]) + cigar_len(right_cigar, CONSUMES_READ) - int(primary[0][POS])
        for num, op in string2cigar(primary[0][CIGAR]):
            if op in "IS":
                size += num
            elif op in "DN":
                size -= num
    
        if size < 0 or (max_fragment_size and size > max_fragment_size):
            size = 0
        
    if size:
        stats["total"][size] += 1
        if contigs:
            rnames = set(segment[RNAME] for segment in read)
            for rname in rnames & contigs:
                stats[rname][size] += 1



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', help="Input sam file, must be sorted by name.")
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="stats_file", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--max-fragment-size", help="Maximum fragment size to be considered a genuine read pair.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-r", "--rnames", help="Reference sequence names over which to calculate fragment size distributions.", default=argparse.SUPPRESS)
    parser.add_argument("-o", "--output", help="Output pdf.", dest="output", default=argparse.SUPPRESS)
    parser.add_argument("-n", "--name", help="Sample name.", default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        do_sizing(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

