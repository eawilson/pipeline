import pdb
import argparse
import json
from collections import defaultdict
from statistics import mean



def viral_copies(statistics, targets, exclude=""):
    
    with open(statistics, "rt") as f_in:
        stats = json.load(f_in)

    targets = set(targets.split())
    exclude = set(exclude.split())
    
    baseline = []
    depths = defaultdict(list)
    for target, depth in stats.get("fragments_per_target", ()):
        target = "_".join(target.split("_")[:-1])
        if target in targets:
            depths[target].append(depth)
        elif target not in exclude:
            baseline.append(depth)
    
    baseline = mean(baseline)
    for target, depth in depths.items():
        stats["viral_copies_per_cell"][target] = mean(depth) / baseline
        
    with open(statistics, "wt") as f_out:
        json.dump(stats, f_out, sort_keys=True, indent=4)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stats', help="Input stats file.", dest="statistics")
    parser.add_argument("-t", "--targets", help="Viral targets.", required=True)
    parser.add_argument("-e", "--exclude", help="Targets to exclude when calculating baseline.", default=argparse.SUPPRESS)
    
    args = parser.parse_args()
    try:
        viral_copies(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

