import pdb
import argparse
import sys
import os
from collections import Counter, defaultdict

#from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
#from matplotlib.figure import Figure
#from matplotlib.patches import Rectangle
#from matplotlib.backends.backend_pdf import PdfPages

from .utils import save_stats


nucleotide = {"A": "R", "G": "R", "C": "Y", "T": "Y"}
REF = 3
ALT = 4
FORMAT = 8
SAMPLE1 = 9

def vcf_stats(vcf_path, stats_file="stats.json", output="vafs.pdf", name=""):
    
    if not name:
        name = os.path.basename(vcf_path).split(".")[0]

    total = 0
    ti = 0
    tv = 0
    snps = 0
    indels = 0
    vafs = []
    with open(vcf_path, "rt") as f_in:
        for variant in f_in:
            if variant.startswith("#") or variant.startswith('"#'):
                continue
            
            variant = variant.split("\t")
            total += 1
            try:
                if nucleotide[variant[REF]] == nucleotide[variant[ALT]]:
                    ti += 1
                else:
                    tv += 1
            except KeyError:
                pass
            
            if len(variant[REF]) != len(variant[ALT]):
                indels += 1
            else:
                snps += 1
            
            #keys = variant[FORMAT].split(":")
            #vals = variant[SAMPLE1].strip().split(":")
            #vaf = dict(zip(keys, vals))["FREQ"]
            #vafs.append(float(vaf.rstrip("%")) / 100)
    
    
    #pdf = PdfPages(output)
    #figure = Figure(figsize=(11.69,8.27))
    #FigureCanvasPdf(figure)
    #ax = figure.gca()
    
    #ax.hist(vafs, bins=100)
    #ax.get_yaxis().set_visible(False)
    #ax.set_xlabel("Variant Allele Frequency", fontsize=10)
    #ax.set_ylabel("Frequency", fontsize=10)
    #ax.set_title(sample, fontsize=12)
    
    #pdf.savefig(figure)
    #pdf.close()
    
    save_stats(stats_file, {"variants": {"total": total,
                                         "snp/indel": float(snps) / (indels or 1),
                                         "ti/tv": float(ti) / (tv or 1)}})



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf_path', help="Input vcf file.")
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="stats_file", default=argparse.SUPPRESS)
    parser.add_argument("-o", "--output", help="Output pdf.", dest="output", default=argparse.SUPPRESS)
    parser.add_argument("-n", "--name", help="Sample name.", default=argparse.SUPPRESS)
    
    args = parser.parse_args()
    try:
        vcf_stats(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

