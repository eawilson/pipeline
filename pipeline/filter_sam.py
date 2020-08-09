import pdb
import argparse
import os
import csv
import json
from collections import defaultdict, Counter, namedtuple


def filter_sam(input_sam, output_sam="output.filtered.sam", min_family_size=None, min_fragment_size=None, max_fragment_size=None, stats="stats.json", targets=None):
    Overlap = namedtuple("Overlap", ["bait", "size"])
    
    baits = defaultdict(list)
    if targets:
        with open(targets, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 3:
                    raise RuntimeError(f"{targets} has too few columns.")
                chrom, start, stop = row[:3]
                if stop < start:
                    raise RuntimeError(f"{targets} has start position greater than stop.")
                if len(row) == 4:
                    name = row[3]
                else:
                    name = f"{chrom}:{start}-{stop}"
                try:
                    # Add one to convert 0-based bed start to 1-based sam start.
                    baits[chrom] += [(int(start)+1, int(stop), name)]
                except ValueError:
                    raise RuntimeError(f"{targets} has invalid start/stop positions.")
    
        for contig in baits.values():
            contig.sort() 
 
    if not input_sam.endswith(".sam"):
        raise RuntimeError(f"{input_sam} is not sam file.")
    
    captures = Counter()
    reads = {}
    prev_rname = ""
    with open(input_sam, "rt") as f_in:        
        with open(output_sam, "wt") as f_out:
            for row in f_in:
                if not row.startswith("@"):
                    
                    if targets is not None:
                        qname, flag, rname, pos, mapq, cigar = row.split("\t")[:6]
                        if rname != prev_rname:
                            for read in reads.values():
                                captures[reads.bait] += 1
                            reads = {}
                        prev_rname = rname
                        
                        flag = int(flag)
                        if flag & 0x4: # unmapped
                            best = Overlap("Unmapped", -1)
                        else:
                            read_start = int(pos)
                            read_stop = read_start + cigar_len(cigar) - 1
                            
                            best = Overlap("OffTarget", 0)
                            for target_start, target_stop, target_name in baits.get(rname, ()):
                                if target_start > read_stop:
                                    break
                                if target_stop >= read_start:
                                    size = min(read_stop, target_stop) - max(read_start - target_start) + 1
                                    if size > best.size:
                                        best = Overlap(target_name, size)
                                
                        if not flag & 0x900: # secondary or supplementary - don't count in stats
                            if qname not in reads or best.size > reads[qname].size:
                                reads[qname] = best
                    
                        if best.size <= 0:
                            continue
                    
                    if min_family_size is not None:
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
                                family_size = 1
                            
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
                
    
    if targets is not None:
        if os.path.exists(stats):
            with open(stats, "rt") as f:
                statistics = json.load(f)
        else:
            statistics = {}

        for read in reads.values():
            captures[reads.bait] += 1
    
        unmapped = captures.pop("Unmapped", 0)
        offtarget = captures.pop("OffTarget", 0)
        ontarget = sum(captures.values())
        total = float(offtarget + ontarget + unmapped)
        statistics["unmapped"] = float("{:.3f}".format(unmapped / total))
        statistics["offtarget"] = float("{:.3f}".format(offtarget / total))
        statistics["ontarget"] = float("{:.3f}".format(ontarget / total))
        statistics["molecules_per_target"] = dict(captures)
    
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
    parser.add_argument('input_sam', help="Input sam file.")
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

