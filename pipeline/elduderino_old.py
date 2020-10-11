import pdb
import argparse
import json
import sys
import csv
from collections import defaultdict, Counter
from itertools import chain
from multiprocessing import Process, Queue
try:
    from contextlib import nullcontext
except ImportError: # python <3.7
    class nullcontext(object):
        def __enter__(self):
            return None

        def __exit__(self, *excinfo):
            pass

QNAME = 0
FLAG = 1
RNAME = 2
POS = 3
MAPQ = 4
CIGAR = 5
TLEN = 8
SEQ = 9
QUAL = 10

UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
RC = 0x10
READ1 = 0x40
READ2 = 0x80
SECONDARY = 0X100
FILTERED = 0x200
SUPPPLEMENTARY = 0x800
NOT_PRIMARY_AND_MAPPED = UNMAPPED | MATE_UNMAPPED | SECONDARY | SUPPPLEMENTARY



def dedupe_process(input_queue, output_queue, stats_queue, dedupe_func, stats, **details):
    while True:
        size_family = input_queue.get()
        if size_family is None:
            break
        output_queue.put(dedupe_func(size_family, stats=stats, **details))
    stats_queue.put(stats)
    


def filewriter_process(output_queue, output_sam):
    with open(output_sam, "wt") as f_out:
        while True:
            text = output_queue.get()
            if text is None:
                break
            f_out.write(text)



def read_bed(targets):
    baits = defaultdict(list)
    with open(targets, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            try:
                chrom = row[0]
                # Add one to convert 0-based bed start to 1-based sam start
                start = int(row[1]) + 1
                stop = int(row[2])
            except (ValueError, IndexError):
                row = " ".join(row)
                sys.exit(f'Unable to parse "{row}" within {targets} bedfile')
            if stop < start:
                sys.exit(f"{targets} has invalid start/stop positions")
            try:
                name = "{}_{}".format(row[3], f"{chrom}:{start}-{stop}")
            except IndexError:
                name = f"{chrom}:{start}-{stop}"
            baits[chrom] += [(start, stop, name)]
    for contig in baits.values():
        contig.sort()
    return baits



def unpack_cigar(cigar_string):
    cigar_tuples = []
    num = ""
    for character in cigar_string:
        if character.isnumeric():
            num += character
        else:
            try:
                cigar_tuples.append((int(num), character))
            except ValueError:
                sys.exit(f"Malformed cigar string {cigar_string}")
            num = ""
    return tuple(cigar_tuples)



def repack_cigar(cigar_tuples):
    return "".join(str(v) for v in chain(*cigar_tuples))



def reference_len(cigar):
    length = 0
    for num, op in cigar:
        if op in "MD=XN":
            length += num
    return length



def elduderino(input_sam, output_sam="output.deduped.sam", statistics="stats.json", umi=None, min_fragment_size=None, max_fragment_size=None, min_family_size=2, targets=None, threads=1):

    if umi == "thruplex":
        dedupe_func = dedupe#_umi_inexact
    elif umi in ("thruplex_hv", "prism"):
        dedupe_func = dedupe_umi_exact
    elif umi is None:
        dedupe_func = dedupe
    else:
        sys.exit(f"'{umi}' is not a known UMI type")
    
    details = {"targets": read_bed(targets) if targets is not None else {},
               "min_fragment_size": min_fragment_size,
               "max_fragment_size": max_fragment_size,
               "min_family_size": min_family_size}
    
    stats = {"fragment_sizes": Counter(),
             "family_sizes": Counter()}
    if targets is not None:
        stats.update({"fragments_per_target": Counter({val[2]: 0 for val in chain(*details["targets"].values())}),
                      "ontarget": 0,
                      "offtarget": 0})
    
    
    multithreaded = (threads > 1)
    with (open(output_sam, "wt") if not multithreaded else nullcontext()) as f_out:
        
        if multithreaded:
            input_queue = Queue()
            output_queue = Queue()
            stats_queue = Queue()
            workers = []
            for i in range(threads - 1):
                workers.append(Process(target=dedupe_process, args=(input_queue, output_queue, stats_queue, dedupe_func, stats), kwargs=details))
            workers.append(Process(target=filewriter_process, args=(output_queue, output_sam)))
            for worker in workers:
                worker.start()
            dedupe_family = input_queue.put
            write = output_queue.put
            
        else:
            def dedupe_family(size_family):
                f_out.write(dedupe_func(size_family, stats=stats, **details))
            write = f_out.write

        
        unpaired = {}
        
        current_rname = ""
        paired = defaultdict(list)
        max_pos = 0
        paired2 = defaultdict(list)
        max_pos2 = 0

        with open(input_sam, "rt") as f_in:
            previous_rnames = set()
            sort_check_rname = ""
            sort_check_pos = 0
            reader = csv.reader(f_in, delimiter="\t")
            for read in reader:
                if read[0].startswith("@"):
                    write("\t".join(read))
                    write("\n")
                    continue
                
                qname = read[QNAME]
                flag = read[FLAG] = int(read[FLAG])
                rname =  read[RNAME]
                pos = read[POS] = int(read[POS])
                
                if rname == sort_check_rname:
                    if pos < sort_check_pos:
                        sys.exit("SAM file must be sorted by position")
                elif sort_check_rname in previous_rnames:
                    sys.exit("SAM file must be sorted by position")
                else:
                    previous_rnames.add(sort_check_rname)
                    sort_check_rname = rname
                sort_check_pos = pos

                if flag & NOT_PRIMARY_AND_MAPPED:
                    continue
                
                if not((flag & READ1) ^ (flag & READ2)):
                    sys.exit("Invalid segment SAM flags")
                
                cigar = read[CIGAR] = unpack_cigar(read[CIGAR])
                
                try:
                    mate = unpaired.pop(qname)
                except KeyError:
                    unpaired[qname] = read
                    continue
                
                read_rc = flag & RC
                read_start = pos if not read_rc else pos + reference_len(cigar) - 1
                mate_rc = mate[FLAG] & RC
                mate_start = mate[POS] if not mate_rc else mate[POS] + reference_len(mate[CIGAR]) - 1
                location = (mate[RNAME], mate_start, mate_rc, rname, read_start, read_rc)
                
                if rname == current_rname and pos <= max_pos:
                    if location in paired:
                        paired[location].append((mate, read))
                    
                    else:
                        paired2[location].append((mate, read))
                        if pos > max_pos2:
                            max_pos2 = pos
                    
                else:
                    if multithreaded and not all(worker.is_alive() for worker in workers):
                        sys.exit()
                    
                    for size_family in paired.values():
                        dedupe_family(size_family)
                    paired = paired2
                    max_pos = max_pos2
                    paired2 = defaultdict(list)
                    max_pos2 = 0
                    
                    paired[location].append((mate, read))
                    if rname != current_rname or pos > max_pos:
                        current_rname = rname
                        max_pos = pos
                    
            for size_family in chain(paired.values(), paired2.values()):
                dedupe_family(size_family)


        if multithreaded:
            for worker in workers[:-1]:
                input_queue.put(None)
            for worker in workers[:-1]:
                _stats = stats_queue.get()
                for key, val in _stats.items():
                    if isinstance(val, Counter):
                        for k, v in val.items():
                            stats[key][k] += v
                    else:
                        stats[key] += _stats[key]
            output_queue.put(None)
            for worker in workers:
                worker.join()

    
    try:
        with open(statistics, "rt") as f:
            old_stats = json.load(f)
    except OSError:
        old_stats = {}
    old_stats.update(stats)
    with open(statistics, "wt") as f:
        json.dump(old_stats, f, sort_keys=True, indent=4)



def dedupe_umi_exact(size_family, **kwargs):
    umi_families = defaultdict(list)
    for pair in size_family:
        for tag in pair[0][11:]:
            if tag.startswith("RX:Z:"):
                umi_families[tag].append(pair)
                break
        else:
            sys.exit("Missing RX tags")
    #if len(size_family) > 1:
        #for umi, family in umi_families.items():
            #print(umi.split(":")[-1], len(family))
        #print("")
    return "".join([dedupe(family, **kwargs) for family in umi_families.values()])



def dedupe_umi_inexact(size_family, **kwargs):
    umi_pairs = []
    #for pair in size_family:
        #for tag in pair[0][11:]:
            #if tag.startswith("RX:Z:"):
                #left, right = tag[5:].split("-")
                #umi_pairs.append([left, right, pair])
                #break
        #else:
            #sys.exit("Missing RX tags")
    
    
    #while umi_pairs:
        #remaining = []
        #left, right, pair = umi_pairs.pop()
        #family = [pair]
        #left_umis = set([left])
        #right_umis = set([right])
        #for pair in umi_pairs:
            #if pair[0] in left_umis or pair[1] in right_umis:
            
            
            
            #else:
                
        
        
        
        
    #umi_families = defaultdict(list)
    #for pair in family:
        #for tag in pair[11:]:
            #if tag.startswith("RX:Z:"):
                #umi_families[tag[5:]].append(pair)
                #break
        #else:
            #sys.exit("Missing RX tags")
    #for sub_family in umi_families.values():
        #dedupe(sub_family)



def dedupe(family, stats, targets, min_family_size, min_fragment_size, max_fragment_size):
    passed = True
    family_size = len(family)
    stats["family_sizes"][family_size] += 1
    
    min_size = (family_size * 6 // 10) + 1
    cigar_families = defaultdict(list)
    for pair in family:
        cigar_families[(pair[0][CIGAR], pair[1][CIGAR])].append(pair)
    family = sorted(cigar_families.values(), key=lambda x:len(x))[-1]
    family_size = len(family)
    if family_size < min_size or family_size < min_family_size:
        passed = False
    
    
    # The segments are the wrong way around, likely due to readthrough.
    # Swap them and trim the readthrough.
    #if left[RNAME] == right[RNAME] and left[FLAG] & RC and not(right[FLAG] & RC):
        #for pair in family:
            #pair[1], pair[0] = pair[0], pair[1]
        
    
    first_pair = family[0]
    left, right = first_pair
    left_rname = left[RNAME]
    right_rname = right[RNAME]
    left_rc = left[FLAG] & RC
    right_rc = right[FLAG] & RC
    left_start = left[POS]
    right_start = right[POS]
    if left_rname == right_rname and not(left_rc) and right_rc:
        #print("PASS", "Read1" if left[FLAG] & READ1 else "Read2", "rvs" if left[FLAG] & RC else "fwd", "Read1" if right[FLAG] & READ1 else "Read2", "rvs" if right[FLAG] & RC else "fwd", left[TLEN], right[TLEN], left[CIGAR], right[CIGAR])
        right_stop = right_start + reference_len(right[CIGAR]) - 1
        fragment_size = right_stop - left_start + 1
        
        if targets:
            best_overlap = 0
            target = None
            for start, stop, name in targets.get(left_rname, ()):
                if start > right_stop:
                    break
                if stop >= left_start:
                    overlap = min(right_stop, stop) - max(left_start, start) + 1
                    if overlap > best_overlap:
                        best_overlap = overlap
                        target = name
            if target is not None:
                stats["fragments_per_target"][target] += 1
                stats["ontarget"] += 1
                passed = True
            else:
                stats["offtarget"] += 1
                passed = False
    
    elif left_rname == right_rname and left_rc and not(right_rc):
        pass
        fragment_size = 0
        if targets:
            passed = False
        
        
        
    
    else:
        #print("{}:{}-{}  {}:{}-{}".format(left_rname, left_start, "<" if left_rc else ">", right_rname, right_start, "<" if right_rc else ">"))
        #print(left[CIGAR], right[CIGAR])
        #print("FAIL", "Read1" if left[FLAG] & READ1 else "Read2", "rvs" if left[FLAG] & RC else "fwd", "Read1" if right[FLAG] & READ1 else "Read2", "rvs" if right[FLAG] & RC else "fwd", left[TLEN], right[TLEN], left[CIGAR], right[CIGAR])
        #print(left[TLEN], right[TLEN], "\n")
        fragment_size = 0
        if targets:
            passed = False
    stats["fragment_sizes"][fragment_size] += 1
    
    if not passed:
        return ""
    
    if fragment_size:
        # We are looking for the position in the left read that overlies the
        # first consumed base in the right read.
        
        # Find the first right read position that consumes reference.
        right_read_pos = 0
        right_ref_pos = right_start
        for num, op in right[CIGAR]:
            if op in "SI":
                right_read_pos += num
            elif op not in "HP":
                break
        
        # Find the left read position that overlies pos for the right read
        # on the reference.
        left_read_pos = 0
        left_ref_pos = left_start
        num = ""
        for num, op in left[CIGAR]:
            if op in "MDN=X":
                if left_ref_pos + num > right_ref_pos:
                    num = right_ref_pos - left_ref_pos
                left_ref_pos += num
            if op in "MIS=X":
                left_read_pos += num
            if left_ref_pos == right_ref_pos:
                break
        
        # The two segments overlap and can threfore be corrected.
        if left_ref_pos == right_ref_pos:
            overlap_len = min(len(left[SEQ]) - left_read_pos, len(right[SEQ]) - right_read_pos)
            for left, right in family:
                left_seq = list(left[SEQ])
                left_qual = list(left[QUAL])
                right_seq = list(right[SEQ])
                right_qual = list(right[QUAL])

                for i in range(overlap_len):
                    if left_seq[left_read_pos+i] != right_seq[right_read_pos+i]:
                        if ord(left_qual[left_read_pos+i]) > ord(right_qual[right_read_pos+i]) + 10:
                            right_seq[right_read_pos+i] = left_seq[left_read_pos+i]
                            right_qual[right_read_pos+i] = left_qual[left_read_pos+i]
                        elif ord(right_qual[right_read_pos+i]) > ord(left_qual[left_read_pos+i]) + 10:
                            left_seq[left_read_pos+i] = right_seq[right_read_pos+i]
                            left_qual[left_read_pos+i] = right_qual[right_read_pos+i]
                        else:
                            right_seq[right_read_pos+i] = "N"
                            right_qual[right_read_pos+i] = "!"
                            left_seq[left_read_pos+i] = "N"
                            left_qual[left_read_pos+i] = "!"

                left[SEQ] = "".join(left_seq)
                left[QUAL] = "".join(left_qual)
                right[SEQ] = "".join(right_seq)
                right[QUAL] = "".join(right_qual)
                
    
    if family_size > 1:
        too_few = family_size * 6 // 10
        for read in range(2):
            #for pair in family:
                #print(pair[read][SEQ])
            #print("")
            
            consensus_seq = []
            consensus_qual = []
            for i in range(len(first_pair[read][SEQ])):
                bases = Counter()
                quals = defaultdict(lambda:"!")
                for pair in family:
                    base = pair[read][SEQ][i]
                    qual = pair[read][QUAL][i]
                    bases[base] += 1
                    if qual > quals[base]:
                        quals[base] = qual
                
                base, count = sorted(Counter(bases).items(), key=lambda x:x[1])[-1]
                if count > too_few:
                    consensus_seq.append(base)
                    consensus_qual.append(quals[base])
                else:
                    consensus_seq.append("N")
                    consensus_qual.append("!")

            first_pair[read][SEQ] = "".join(consensus_seq)
            first_pair[read][QUAL] = "".join(consensus_qual)
            first_pair[read][MAPQ] = str(max(int(pair[read][MAPQ]) for pair in family))
            #print(first_pair[read][SEQ], "\n")
    
    for read in range(2):
        first_pair[read][POS] = str(first_pair[read][POS])
        first_pair[read][FLAG] = str(first_pair[read][FLAG])
        first_pair[read][CIGAR] = repack_cigar(first_pair[read][CIGAR])
    
    return "{}\n{}\n".format("\t".join(first_pair[0]), "\t".join(first_pair[1]))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', help="Input sam file.")
    parser.add_argument("-o", "--output", help="Output sam file.", dest="output_sam", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--min-family-size", help="Minimum family size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-f", "--min-fragment-size", help="Minimum fragment size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-F", "--max-fragment-size", help="Maximum fragment size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="statistics", default=argparse.SUPPRESS)
    parser.add_argument("-b", "--bed", help="Bed file of on-target regions.", dest="targets", default=argparse.SUPPRESS)
    parser.add_argument("-u", "--umi", help="UMI type, allowed = thruplex, thruplex_hv, prism.", default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        elduderino(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

