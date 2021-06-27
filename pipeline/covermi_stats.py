import pdb
import argparse
import sys
import os
from collections import Counter, defaultdict
from covermi import Panel, Cov, Plot

from .utils import save_stats


def covermi_stats(bam_path, panel_path, output_path="coverage.pdf", stats_file="stats.json", name=""):
    if not name:
        name = os.path.basename(bam_path).split(".")[0]

    panel = Panel(panel_path)
    if "targets" in panel:
        roi = panel.targets
    elif "exons" in panel:
        roi = panel.exons
    else:
        return
    cov = Cov(bam_path)
    Plot(coverage=cov, panel=panel, depth=None, title=name, output=output_path)
    
    stats = {"coverage": {},
             "coverage_by_gene": defaultdict(dict)}
    for depth in (30, 100, 500, 1000, 2000):
        for i in cov.calculate(roi, depth):
            stats["coverage_by_gene"][f"{depth}x"][i.name] = int(i.percent_covered)
            stats["coverage_by_gene"]["mean_depth"][i.name] = int(i.depth)
    
        i = cov.calculate(roi, depth, name="Total")
        stats["coverage"][f"{depth}x"] = int(i.percent_covered)
        stats["coverage"]["mean_depth"] = int(i.depth)
    
    save_stats(stats_file, stats)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_path', help="Input bam file.")
    parser.add_argument("-p", "--panel", help="Covermi panel.", dest="panel_path", required=True)
    parser.add_argument("-o", "--output", help="Output coverage plot.", dest="output_path", default=argparse.SUPPRESS)
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="stats_file", default=argparse.SUPPRESS)
    parser.add_argument("-n", "--name", help="Sample name.", default=argparse.SUPPRESS)
    
    args = parser.parse_args()
    try:
        covermi_stats(**vars(args))
    except OSError as e:
        raise
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

