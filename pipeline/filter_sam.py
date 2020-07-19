import pdb
import argparse
import os
import csv
import json
from collections import defaultdict, Counter, namedtuple


def filter_sam(input_sam, output_sam="output.filtered.sam", min_family_size=None, min_fragment_size=None, max_fragment_size=None, stats="stats.json", targets=None):
    Cluster = namedtuple("Cluster", ["target_name", "overlap", "family_size"])
    
    input_sam = input_sam[0]
    if not input_sam.endswith(".sam"):
        raise RuntimeError(f"{input_sam} is not sam file.")
    
    contigs = defaultdict(list)
    if targets is not None:
        if not os.path.exists(targets):
            raise RuntimeError(f"{targets} does not exist.")
    
        with open(targets, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 3:
                    raise RuntimeError(f"{targets} has too few columns.")
                chrom, start, stop = row[:3]
                if len(row) == 4:
                    name = row
                else:
                    name = f"{chrom}:{start}-{stop}"
                try:
                    contigs[chrom] += [(int(start), int(stop), name)]
                except ValueError:
                    raise RuntimeError(f"{targets} has invalid start/stop positions.")
        
        for val in contigs.values():
            val.sort()
    
    clusters = {}
    with open(input_sam, "rt") as f_in:        
        with open(output_sam, "wt") as f_out:
            for row in f_in:
                
                if not row.startswith("@"):
                    
                    start = row.find("\tYF:i:")
                    if start == -1:
                        family_size = 1
                    else:
                        stop = row.find("\t", start + 6)
                        if stop == -1:
                            stop = len(row)
                        try:
                            family_size = int(row[start + 6:stop])
                        except ValueError:
                            family_size = 1 # ?print warning
                    
                    if targets is not None:
                        qname, flag, rname, pos, mapq, cigar = row.split("\t")[:6]
                        flag = int(flag)
                        if flag & 0x900 or rname not in contigs: # secondary or supplementary
                            continue
                        if flag & 0x4: # unmapped
                            if qname not in clusters:
                                clusters[qname] = Cluster(None, None, family_size)
                            continue
                        
                        read_start = int(pos) - 1
                        read_stop = read_start + cigar_len(cigar)
                        
                        best_overlap = 0
                        best_name = None
                        for target_start, target_stop, target_name in contigs[rname]:
                            if target_start >= read_stop:
                                break
                            if target_stop > read_start:
                                if target_start < read_start:
                                    overlap = (read_stop if read_stop < target_stop else target_stop) - (read_start if read_start > target_start else target_start)
                                    if overlap > best_overlap:
                                        best_overlap = overlap
                                        best_name = target_name
                                    elif overlap < 1:
                                        raise RuntimeError("Logic error calculating overlap.")
                        if best_name is None:
                            if qname not in clusters: 
                                clusters[qname] = Cluster(None, None, family_size)
                            continue
                        else:
                            if qname in clusters:
                                cluster = clusters[qname]
                                if cluster.target_name is None:
                                    clusters[qname] = Cluster(best_name, best_overlap, family_size)
                                elif best_overlap > cluster.overlap:
                                     clusters[qname] = Cluster(best_name, best_overlap, family_size)
                            else:
                                clusters[qname] = Cluster(best_name, best_overlap, family_size)

                            
                    if min_family_size is not None:
                        if family_size < min_family_size:
                            continue

                    
                    if min_fragment_size is not None or max_fragment_size is not None:
                        start = row.find("\tYS:i:")
                        if start == -1:
                            continue
                        stop = row.find("\t", start + 6)
                        if stop == -1:
                            stop = len(row)                            
                        try:
                            fragment_size = int(row[start + 6:stop])
                        except ValueError:
                            continue
                        if (min_fragment_size is not None and fragment_size < min_fragment_size) or (max_fragment_size is not None and fragment_size > max_fragment_size):
                            continue
                    
                f_out.write(row)
                
    
    if targets is not None and clusters:
        if os.path.exists(stats):
            with open(stats, "rt") as f:
                statistics = json.load(f)
        else:
            statistics = {}

        dest = Counter()
        offtarget = 0
        ontarget = 0
        for cluster in clusters.values():
            if cluster.target_name is None:
                offtarget += cluster.family_size
            else:
                dest[cluster.target_name] += 1
                ontarget += cluster.family_size
    
        statistics["offtarget"] = float("{:.3f}".format(offtarget / float(offtarget + ontarget)))
        if dest:
            statistics["molecules_per_target"] = dict(dest)
    
        with open(stats, "wt") as f:
            json.dump(statistics, f, sort_keys=True, indent=4)
    


def cigar_len(cigar):
    num = ""
    length = 0
    for character in cigar:
        if character.isnumeric():
            num += character
        else:
            if character in "MD=XN":
                length += int(num)
            elif character not in "ISHP":
                raise RuntimeError("Malformed cigar string.")
            num = ""
    return length



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', nargs=1, help="Input sam file.")
    parser.add_argument("-o", "--output", help="Output sam file.", dest="output_sam", default=argparse.SUPPRESS)
    parser.add_argument("-f", "--min-family-size", help="Minimum family size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-z", "--min-fragment-size", help="Minimum fragment size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-Z", "--max-fragment-size", help="Maximum fragment size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-s", "--stats", help="Statistics file.", default=argparse.SUPPRESS)
    parser.add_argument("-t", "--targets", help="Bed file of on-target regions.", default=argparse.SUPPRESS)
    args = parser.parse_args()
    filter_sam(**vars(args))



if __name__ == "__main__":
    main()

