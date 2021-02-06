import pdb
import argparse
import sys
import os
import glob
from collections import Counter, defaultdict
from multiprocessing import Process, Queue
from collections.abc import Mapping

from covermi import bed, Gr, Entry

from .utils import run, save_stats, string2cigar

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
SEC_OR_SUP = SECONDARY | SUPPPLEMENTARY
BOTH_UNMAPPED = UNMAPPED | MATE_UNMAPPED

CONSUMES_REF = "MDN=X"
CONSUMES_READ = "MIS=X"

LEFT = 0
RIGHT = 1

RCOMPLEMENT = str.maketrans("ATGC", "TACG")



def mean(nums):
    return (sum(nums) * 1.0) / len(nums)



def worker_process(input_queue, output_queue, stats_queue, stats, **details):
    while True:
        read = input_queue.get()
        if read is None:
            break
        output = _filter_read(read, stats=stats, **details)
        if output:
            output_queue.put(output)
    stats_queue.put(stats)
    


def filewriter_process(output_queue, output_file):
    with open(output_file, "wt") as f_out:
        while True:
            output = output_queue.get()
            if output is None:
                break
            f_out.write(output)



def cigar_len(cig, ops):
    return sum(num for num, op in cig if op in ops)



def filter_sam(input_sam,
               output_file="output.filtered.sam",
               stats_file="stats.json",
               targets="",
               max_fragment_size=1000,
               filter_fragments_shorter_than=0,
               filter_fragments_longer_than=0,
               #filter_mapq_less_than=0,
               #filter_low_qual_bases_greater_than=0,
               fragment_sizes_per_target="",
               threads=0):

    threads = 1
    if not threads:
        threads = int(run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip())

    if targets and os.path.isdir(targets):
        beds = glob.glob(f"{targets}/*.bed")
        if len(beds) == 0:
            sys.exit(f'"{targets}" does not contain a bedfile')
        elif len(beds) > 1:
            sys.exit(f'"{targets}" contains multiple bedfiles')
        targets =  beds[0]
    
    targets = Gr(bed(targets or ()))
    for target in targets:
        loc = f"{target.chrom}:{target.start}-{target.stop}"
        target.name = loc if target.name == "." else f"{target.name}_{loc}"
    details = {"targets": targets,
               "max_fragment_size": max_fragment_size,
               "filter_fragments_shorter_than": filter_fragments_shorter_than,
               "filter_fragments_longer_than": filter_fragments_shorter_than}
    
    stats = {"fragment_sizes": Counter()}
    if targets:
        stats.update({"fragments_per_target": {bait.name: 0 for bait in details["targets"]},
                      "fragment_sizes_per_target": defaultdict(Counter),
                      "ontarget": 0,
                      "offtarget": 0})
    
    
    multithreaded = (threads > 1)
    with (open(output_file, "wt") if not multithreaded else nullcontext()) as f_out:
        
        if multithreaded:
            input_queue = Queue()
            output_queue = Queue()
            stats_queue = Queue()
            workers = []
            for i in range(threads - 1):
                workers.append(Process(target=worker_process, args=(input_queue, output_queue, stats_queue, stats), kwargs=details))
            workers.append(Process(target=filewriter_process, args=(output_queue, output_file)))
            for worker in workers:
                worker.start()
            filter_read = input_queue.put
            write = output_queue.put
            
        else:
            def filter_read(read):
                output = _filter_read(read, stats=stats, **details)
                if output:
                    f_out.write(output)
            write = f_out.write

        with open(input_sam, "rt") as f_in:
            current_qname = ""
            read = []
            for row in f_in:
                if row.startswith("@"):
                    write(row)
                    continue
                
                segment = row.split("\t")
                qname = segment[QNAME]                    
                if qname != current_qname:
                    if multithreaded and not all(worker.is_alive() for worker in workers):
                        sys.exit("Worker thread unexpectedly terminated")
                    filter_read(read)
                    current_qname = qname
                    read = []
                read.append(segment)
                
            filter_read(read)


    if multithreaded:
        for worker in workers[:-1]:
            input_queue.put(None)
        for worker in workers[:-1]:
            _stats = stats_queue.get()
            for key, val in _stats.items():
                if isinstance(val, Mapping):
                    for k, v in val.items():
                        stats[key][k] += v
                else:
                    stats[key] += _stats[key]
        output_queue.put(None)
        for worker in workers:
            worker.join()


    include = set(fragment_sizes_per_target.split())
    sizes = defaultdict(Counter)
    for name, depth in stats["fragment_sizes_per_target"].items():
        for target in set(name[name.find("_")+1:].strip().split(";")):
            target = target.split()[0]
            if target in include:
                sizes[target][depth] += 1
    stats["fragment_sizes_per_target"] = sizes
    
    
    baseline = mean(baseline)
    stats = {"copies_per_cell": {t: mean(d) / baseline for t, d in depths.items()}}
    ontarget = stats.pop("ontarget")
    offtarget = stats.pop("offtarget")
    stats["offtarget"] = float(offtarget) / (ontarget + offtarget)    
    save_stats(stats_file, stats)



def _filter_read(read, stats, targets, max_fragment_size, filter_fragments_shorter_than, filter_fragments_longer_than):
    if not read:
        return ""
    
    primary = []
    non_primary = []
    for segment in read:
        segment[FLAG] = int(segment[FLAG])
        if not segment[FLAG] & SEC_OR_SUP: # Primary read
            primary.append(segment)
        else:
            non_primary.append(segment)

    if len(primary) != 2:
        sys.exit("Multiple primary reads in SAM file")
    
    # Both mapped to same reference and pointing in opposite directions
    # therefore may be a concordant read pair.
    if not (primary[0][FLAG] & UNMAPPED) and not (primary[1][FLAG] & UNMAPPED) and \
        primary[0][RNAME] == primary[1][RNAME] and \
        primary[0][FLAG] & RC != primary[1][FLAG] & RC:
        
        if primary[0][FLAG] & RC:
            primary = primary[::-1]
        
        right_cigar = string2cigar(primary[1][CIGAR])
        start = int(primary[0][POS])
        stop = int(primary[1][POS]) + cigar_len(right_cigar, CONSUMES_REF) - 1
        size = (int(primary[1][POS]) + cigar_len(right_cigar, CONSUMES_READ)) - start - 1
        for num, op in string2cigar(primary[0][CIGAR]):
            if op in "IS":
                size += num
            elif op in "DN":
                size -= num
        if size < 0 or (max_fragment_size and size > max_fragment_size):
            size = 0
        
    else:
        size = 0
    
    stats["fragment_sizes"][size] += 1
    
    if targets:
        if size:
            segments = [Entry(primary[0][RNAME], start, stop)]
            remaining = non_primary
        else:
            segments = []
            remaining = read
            
        for segment in remaining:
            start = int(segment[POS])
            segments.append(Entry(segment[RNAME], start, start + cigar_len(string2cigar(segment[CIGAR]), CONSUMES_REF)))
        
        best_offset = None
        match = None
        for segment in segments:
            for target in targets.touched_by(segment):
                offset = abs(segment.start + segment.stop - target.start - target.stop)
                try:
                    if offset > best_offset:
                        continue
                except TypeError:
                    pass
                best_offset = offset
                match = target.name
        
        if match:
            stats["fragments_per_target"][match] += 1
            stats["ontarget"] += 1
            if size:
                stats["fragment_sizes_per_target"][match][size] += 1
        else:
            stats["offtarget"] += 1
            return ""
    
    if filter_fragments_shorter_than and (not size or size < filter_fragments_shorter_than):
        return ""
    if filter_fragments_longer_than and (not size or size > filter_fragments_longer_than):
        return ""
    
    for segment in read:
        segment[FLAG] = str(segment[FLAG])
    return "".join("\t".join(segment) for segment in read)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', help="Input sam file, must be sorted by name.")
    parser.add_argument("-o", "--output", help="Output file.", dest="output_file", default=argparse.SUPPRESS)
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="stats_file", default=argparse.SUPPRESS)
    parser.add_argument("-b", "--bed", help="Bed file of on-target regions.", dest="targets", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--max-fragment-size", help="Maximum fragment size to be considered a genuine read pair.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-f", "--filter-fragments-shorter-than", help="Filter fragments shorter than.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-F", "--filter-fragments-longer-than", help="Filter fragments longer than.", type=int, default=argparse.SUPPRESS)
    #parser.add_argument("-M", "--filter-mapq-less-than", help="Filter fragments with mapq less than.", type=int, default=argparse.SUPPRESS)
    #parser.add_argument("-q", "--filter-low-qual-bases-greater-than", help="Filter fragments with greater than low quality bases.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-p", "--fragment-sizes-targets", help="Target names over which to calculate fragment size distribution.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        filter_sam(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

