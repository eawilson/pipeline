import pdb
import argparse
import sys
import json
from collections import Counter, defaultdict
from covermi import Panel, Cov, Plot



def covermi_stats(bam_path, panel_path, output_path="coverage.pdf", statistics="stats.json", sample=""):
    
    panel = Panel(panel_path)
    if "targets" in panel:
        roi = panel.targets
    elif "exons" in panel:
        roi = panel.exons
    else:
        return
    cov = Cov(bam_path)
    Plot(coverage=cov, panel=panel, depth=None, title=sample, output=output_path)
    
    stats = {"coverage": {},
             "coverage_by_gene": defaultdict(dict)}
    for depth in (30, 100, 1000):
        for i in cov.calculate(roi, depth):
            stats["coverage_by_gene"][f"{depth}x"][i.name] = "{:.0f}".format(i.percent_covered)
            stats["coverage_by_gene"]["mean_depth"][i.name] = "{:.0f}".format(i.depth)
    
        i = cov.calculate(roi, depth, name="Total")
        stats["coverage"][f"{depth}x"] = "{:.0f}".format(i.percent_covered)
        stats["coverage"]["mean_depth"] = "{:.0f}".format(i.depth)
    
    try:
        with open(statistics, "rt") as f:
            old_stats = json.load(f)
    except OSError:
        old_stats = {}
    old_stats.update(stats)
    with open(statistics, "wt") as f:
        json.dump(old_stats, f, sort_keys=True, indent=4)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_path', help="Input bam file.")
    parser.add_argument("-p", "--panel", help="Covermi panel.", dest="panel_path", required=True)
    parser.add_argument("-o", "--output", help="Output coverage plot.", dest="output_path", default=argparse.SUPPRESS)
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="statistics", default=argparse.SUPPRESS)
    parser.add_argument("-n", "--sample", help="Sample name.", default=argparse.SUPPRESS)
    
    args = parser.parse_args()
    try:
        covermi_stats(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

