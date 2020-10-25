import pdb
import argparse
import json
from collections import defaultdict
from statistics import mean



def panel_copy_number(statistics, targets, exclude=""):
    
    with open(statistics, "rt") as f_in:
        stats = json.load(f_in)

    targets = set(targets.split())
    exclude = set(exclude.split())
    
    baseline = []
    depths = defaultdict(list)
    for target, depth in stats.get("fragments_per_target", {}).items():
        target = "_".join(target.split("_")[:-1])
        if target in targets:
            depths[target].append(depth)
        elif target not in exclude:
            baseline.append(depth)
    
    baseline = mean(baseline)
    stats["copies_per_cell"] = {}
    for target, depth in depths.items():
        stats["copies_per_cell"][target] = mean(depth) / baseline
        
    with open(statistics, "wt") as f_out:
        json.dump(stats, f_out, sort_keys=True, indent=4)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('statistics', help="Input stats file.")
    parser.add_argument("-t", "--targets", help="Targets over which to calculate copy numbers.", required=True)
    parser.add_argument("-e", "--exclude", help="Targets to exclude from baseline.", default=argparse.SUPPRESS)
    
    args = parser.parse_args()
    try:
        panel_copy_number(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

