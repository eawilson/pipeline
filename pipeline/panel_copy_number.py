import pdb
import argparse
import json
from collections import defaultdict
from statistics import mean

from .utils import save_stats


def panel_copy_number(stats_file, targets):

    with open(stats_file, "rt") as f:
        try:
            fragments_per_target = json.load(f)["fragments_per_target"]
        except (json.JSONDecodeError, KeyError):
            sys.exit(f"Invalid stats file {stats_file}")

    include = set()
    exclude = set()
    for target in targets.split():
        if target.startswith("-"):
            exclude.add(target[1:])
        else:
            include.add(target)
    
    baseline = []
    depths = defaultdict(list)
    for name, depth in fragments_per_target.items():
        for target in set(name[name.find("_")+1:].strip().split(";")):
            target = target.split()[0]
            if target in include:
                depths[target].append(depth)
            elif target not in exclude:
                baseline.append(depth)
    
    baseline = mean(baseline)
    stats = {"copies_per_cell": {t: mean(d) / baseline for t, d in depths.items()}}
    save_stats(stats_file, stats)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stats_file', help="Input stats file.")
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

