import pdb
import argparse
import json
import sys
import csv
from collections import defaultdict, Counter
from itertools import chain
from multiprocessing import Process, Queue
from collections.abc import Sequence

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
RNEXT = 6
PNEXT = 7
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

CONSUMES_REF = "MDN=X"
CONSUMES_READ = "MIS=X"



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



def bed(path):
    baits = defaultdict(list)
    if path:
        with open(path, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                try:
                    chrom = row[0]
                    # Add one to convert 0-based bed start to 1-based sam start
                    start = int(row[1]) + 1
                    stop = int(row[2])
                except (ValueError, IndexError):
                    row = " ".join(row)
                    sys.exit(f'Unable to parse "{row}" within {path} bedfile')
                if stop < start:
                    sys.exit(f"{path} has invalid start/stop positions")
                try:
                    name = "{}_{}".format(row[3], f"{chrom}:{start}-{stop}")
                except IndexError:
                    name = f"{chrom}:{start}-{stop}"
                baits[chrom] += [(start, stop, name)]
        for contig in baits.values():
            contig.sort()
    return baits


# add binary search +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def is_ontarget(targets, *coordinates, stats=None):
    best_offset = None
    match = None
    for rname, start, stop in coordinates:
        for bstart, bstop, bname in targets.get(rname, ()):
            if bstart > stop:
                break
            if bstop >= start:
                offset = abs(start + stop - bstart - bstop)
                try:
                    if offset > best_offset:
                        continue
                except TypeError:
                    pass
                best_offset = offset
                match = bname
    
    if match is not None:
        stats["fragments_per_target"][match] += 1
        stats["ontarget_deduplicated_reads"] += 1
        return True
    else:
        stats["offtarget_deduplicated_reads"] += 1
        return False



class Cigar(Sequence):
    def __init__(self, cigar_string):
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
        if num:
            sys.exit(f"Malformed cigar string {cigar_string}")
        self._cigar = tuple(cigar_tuples)
        self.reference_len = sum(num for num, op in cigar_tuples if op in CONSUMES_REF)

    def __repr__(self):
        return "{}({})".format(type(self).__name__, repr(str(self)))

    def __str__(self):
        return "".join(str(v) for v in chain(*self._cigar))
    
    def __getitem__(self, index):
        return self._cigar[index]
    
    def __len__(self):
        return len(self._cigar)
    
    def __hash__(self):
        return hash(self._cigar)



def elduderino(input_sam,
               output_sam="output.deduped.sam",
               statistics="stats.json",
               umi=None,
               filter_fragments_shorter_than=None,
               filter_fragments_longer_than=None,
               max_fragment_size=None,
               min_family_size=1,
               targets=None,
               threads=1,
               dont_dedupe=False):
    
    if umi == "thruplex":
        dedupe_func = dedupe#_umi_inexact
    elif umi in ("thruplex_hv", "prism"):
        dedupe_func = dedupe_umi_exact
    elif umi is None:
        dedupe_func = dedupe
    else:
        sys.exit(f"'{umi}' is not a known UMI type")
    
    details = {"targets": bed(targets),
               "max_fragment_size": max_fragment_size,
               "filter_fragments_shorter_than": filter_fragments_shorter_than,
               "filter_fragments_longer_than": filter_fragments_longer_than,
               "min_family_size": 1 if dont_dedupe else min_family_size}
    
    stats = {"fragment_sizes": Counter(),
             "family_sizes": Counter()}
    if targets is not None:
        stats.update({"fragments_per_target": Counter({val[2]: 0 for val in chain(*details["targets"].values())}),
                      "ontarget_deduplicated_reads": 0,
                      "offtarget_deduplicated_reads": 0})
    
    
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
                rname = read[RNAME]
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

                if flag & NOT_PRIMARY_AND_MAPPED: #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    continue
                
                if not((flag & READ1) ^ (flag & READ2)):
                    sys.exit("Invalid segment SAM flags")
                
                cigar = read[CIGAR] = Cigar(read[CIGAR])
                
                try:
                    mate = unpaired.pop(qname)
                except KeyError:
                    unpaired[qname] = read
                    continue
                
                if dont_dedupe:
                    dedupe_family([(mate, read)])
                    continue
                
                read_rc = flag & RC
                read_begin = pos if not read_rc else pos + cigar.reference_len - 1
                mate_rc = mate[FLAG] & RC
                mate_begin = mate[POS] if not mate_rc else mate[POS] + mate[CIGAR].reference_len - 1
                location = (mate[RNAME], mate_begin, mate_rc, rname, read_begin, read_rc)
                
                if mate_begin > read_begin:
                    read_begin = mate_begin
                if rname == current_rname and pos <= max_pos:
                    if location in paired:
                        paired[location].append((mate, read))
                    
                    else:
                        paired2[location].append((mate, read))
                        if read_begin > max_pos2:
                            max_pos2 = read_begin
                    
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
                    if rname != current_rname:
                        current_rname = rname
                        max_pos = read_begin
                    elif read_begin > max_pos:
                        max_pos = read_begin
                    
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
    if len(size_family) > 1:
        umi_families = defaultdict(list)
        for pair in size_family:
            for tag in pair[0][11:]:
                if tag.startswith("RX:Z:"):
                    umi_families[tag].append(pair)
                    break
            else:
                sys.exit("Missing RX tags")
        return "".join([dedupe(family, **kwargs) for family in umi_families.values()])
    
    else:
        return dedupe(size_family, **kwargs)


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



def hard_clip(cigar_op):
    return cigar_op[0] if cigar_op[1] == "H" else 0



def dedupe(family, stats, targets, min_family_size, max_fragment_size, filter_fragments_shorter_than, filter_fragments_longer_than):
    passed = True
    family_size = len(family)
    stats["family_sizes"][family_size] += 1
    
    if family_size > 1:
        sixty_percent = (family_size * 6 // 10) + 1
        cigar_families = defaultdict(list)
        for pair in family:
            cigar_families[(pair[0][CIGAR], pair[1][CIGAR])].append(pair)
        family = sorted(cigar_families.values(), key=lambda x:len(x))[-1]
        family_size = len(family)
        if family_size < sixty_percent:
            passed = False
    if family_size < min_family_size and passed:
        passed = False

    readthrough = False
    segments_overlap = False
    fragment_size = 0
    
    first_pair = family[0]
    left, right = first_pair
    left_rname = left[RNAME]
    right_rname = right[RNAME]
    left_rc = left[FLAG] & RC
    right_rc = right[FLAG] & RC
    left_start = left[POS]
    right_start = right[POS]
    left_read_pos = -1
    right_read_pos = 0
    if left_rname == right_rname:
        
        # Do the two segments overlap? Find left_read_pos that overlaps right_start
        ref_pos = left_start - 1
        for num, op in left[CIGAR]:
            if op in CONSUMES_REF:
                if ref_pos + num > right_start:
                    num = right_start - ref_pos
                ref_pos += num
            if op in CONSUMES_READ:
                left_read_pos += num
            if ref_pos == right_start:
                # Subtract non ref consuming bases at the start of right
                # left_read_pos will now overlap the first baes of right
                for num, op in right[CIGAR]:
                    if op in CONSUMES_REF:
                        break
                    if op in CONSUMES_READ:
                        right_read_pos += num
                segments_overlap = True
                break
        
        
        if segments_overlap:
            if left_rc and not right_rc:
                # Inverted read directions therefore readthrough into opposite umi
                readthrough = True
                fragment_size = len(left[SEQ]) - left_read_pos
                if targets and not is_ontarget(targets, (left_rname, right_start, left_start + left[CIGAR].reference_len - 1), stats=stats):
                    passed = False
                
            elif not left_rc and right_rc:
                # Correct read directions therefore no readthrough
                fragment_size = len(right[SEQ]) + left_read_pos
                if targets and not is_ontarget(targets, (left_rname, left_start, right_start + right[CIGAR].reference_len - 1), stats=stats):
                    passed = False

        elif not left_rc and right_rc:
            # Segments don't overlap but are correctly orientated
            fragment_size = right_start - left_start + len(right[SEQ])
            for num, op in left[CIGAR]:
                if op in "IS":
                    fragment_size += num
                elif op in "DN":
                    fragment_size -= num
            if max_fragment_size is not None and fragment_size > max_fragment_size:
                fragment_size = 0
            elif targets and not is_ontarget(targets, (left_rname, left_start, right_start + right[CIGAR].reference_len - 1), stats=stats):
                passed = False
    
        
    stats["fragment_sizes"][fragment_size] += 1
    if not fragment_size and targets:
        coordinates = ((left_rname, left_start, left_start + left[CIGAR].reference_len - 1), (right_rname, right_start, right_start + right[CIGAR].reference_len - 1))
        if not is_ontarget(targets, *coordinates, stats=stats):
            passed = False
    
    
    if filter_fragments_shorter_than is not None and fragment_size >= filter_fragments_shorter_than:
        passed = False
    if filter_fragments_longer_than is not None and fragment_size < filter_fragments_longer_than:
        passed = False

    if not passed:
        return ""
    
    return ""#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    if segments_overlap:
        overlap_len = len(left[SEQ]) - left_read_pos
        for left, right in family:
            left_seq = list(left[SEQ])
            left_qual = list(left[QUAL])
            right_seq = list(right[SEQ])
            right_qual = list(right[QUAL])

            for i in range(overlap_len):
                if left_seq[left_read_pos+i] != right_seq[i]:
                    if ord(left_qual[left_read_pos+i]) > ord(right_qual[i]) + 10:
                        right_seq[i] = left_seq[left_read_pos+i]
                        right_qual[i] = left_qual[left_read_pos+i]
                    elif ord(right_qual[i]) > ord(left_qual[left_read_pos+i]) + 10:
                        left_seq[left_read_pos+i] = right_seq[i]
                        left_qual[left_read_pos+i] = right_qual[i]
                    else:
                        right_seq[i] = "N"
                        right_qual[i] = "!"
                        left_seq[left_read_pos+i] = "N"
                        left_qual[left_read_pos+i] = "!"

            left[SEQ] = "".join(left_seq)
            left[QUAL] = "".join(left_qual)
            right[SEQ] = "".join(right_seq)
            right[QUAL] = "".join(right_qual)
    
    
    if family_size > 1:
        sixty_percent = family_size * 6 // 10
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
                
                base, count = sorted(bases.items(), key=lambda x:x[1])[-1]
                if count > sixty_percent:
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
        first_pair[read][FLAG] = str(first_pair[read][FLAG])
        first_pair[read][POS] = str(first_pair[read][POS])
        first_pair[read][CIGAR] = str(first_pair[read][CIGAR])
    
    return "{}\n{}\n".format("\t".join(first_pair[0]), "\t".join(first_pair[1]))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', help="Input sam file.")
    parser.add_argument("-o", "--output", help="Output sam file.", dest="output_sam", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--min-family-size", help="Minimum family size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-M", "--max-fragment-size", help="Maximum template legth to be considered a genuine read pair.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-f", "--filter-fragments-shorter-than", help="Filter fragments shorter than.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-F", "--filter-fragments-longer-than", help="Filter fragments longer than.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="statistics", default=argparse.SUPPRESS)
    parser.add_argument("-b", "--bed", help="Bed file of on-target regions.", dest="targets", default=argparse.SUPPRESS)
    parser.add_argument("-u", "--umi", help="UMI type, allowed = thruplex, thruplex_hv, prism.", default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-d", "--dont-dedupe", help="Assume sam is already deduped.", action="store_const", const=True, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        elduderino(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

