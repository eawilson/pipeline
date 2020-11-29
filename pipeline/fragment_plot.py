from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages

import pdb
import argparse
import json
import sys
import os


def fragment_plot(stats, output="fragment.sizes.pdf", sample=""):
    
    if not sample:
        sample = os.path.basename(stats).split(".")[0]
    
    with open(stats, "rt") as f:
        try:
            sizes = json.load(f)["fragment_sizes"]
            sizes.pop("0", None)
            size, count = zip(*sorted([(int(x), int(y)) for x, y in sizes.items()]))
        except (json.JSONDecodeError, KeyError, ValueError):
            sys.exit(f"Invalid stats file {stats}")
    
    target = sum(count) // 2
    for median, n in zip(size, count):
        target -= n
        if target <= 0:
            break
    
    pdf = PdfPages(output)
    figure = Figure(figsize=(11.69,8.27))
    FigureCanvasPdf(figure)
    ax = figure.gca()
    
    ax.plot(size, count, "-", color="dodgerblue", linewidth=1)
    ax.get_yaxis().set_visible(False)
    ax.axvline(median, color="black", linewidth=0.5, linestyle=":")
    ax.set_xlabel("Fragment size (bp)", fontsize=10)
    ax.set_title(sample, fontsize=12)
    
    pdf.savefig(figure)
    pdf.close()



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stats', help="Statistics file.")
    parser.add_argument("-o", "--output", help="Output pdf.", dest="output", default=argparse.SUPPRESS)
    parser.add_argument("-n", "--sample", help="Sample name.", default=argparse.SUPPRESS)
    
    args = parser.parse_args()
    try:
        fragment_plot(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()



