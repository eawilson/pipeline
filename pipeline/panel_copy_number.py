import pdb
import argparse
import json
from collections import defaultdict
from statistics import mean



def panel_copy_number(statistics, targets):
    
    with open(statistics, "rt") as f_in:
        stats = json.load(f_in)

    include = set()
    exclude = set()
    for target in targets.split():
        if target.startswith("-"):
            exclude.add(target[1:])
        else:
            include.add(target)
    
    baseline = []
    depths = defaultdict(list)
    for target, depth in stats.get("fragments_per_target", {}).items():
        target = "_".join(target.split("_")[:-1])
        if target in include:
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
    parser.add_argument("-t", "--targets", help="Target names over which to calculate copy numbers. " \
                                                "names preceeded by - will be excluded from baseline.", required=True)
    
    args = parser.parse_args()
    try:
        panel_copy_number(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

