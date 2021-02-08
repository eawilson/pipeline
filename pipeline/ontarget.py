import pdb
import argparse
import sys
import os
import glob
import re
from itertools import chain
from collections import Counter, defaultdict
from multiprocessing import Process, Queue
from collections.abc import Mapping
from statistics import mean

from covermi import bed, Gr, Entry

from .utils import run, save_stats, string2cigar, CONSUMES_REF, CONSUMES_READ

try:
    from contextlib import nullcontext
except ImportError: # python <3.7
    from .utils import nullcontext


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

LEFT = 0
RIGHT = 1



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



def ontarget(input_sam,
             bed_file,
             output_file="output.filtered.sam",
             stats_file="stats.json",
             max_fragment_size=1000,
             retain_offtarget=False,
             cnv = "",
             threads=0):

    threads = 1
    if not threads:
        threads = int(run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip())

    bed_file = (glob.glob(f"{bed_file}/*.bed") + [bed_file])
    if len(bed_file) > 2:
        sys.exit(f'"{targets}" contains multiple bedfiles')
    
    sep = re.compile("[,;]")
    bait2genes= {}
    targets = Gr(bed(bed_file[0]))
    for target in targets:
        genes = set(gene.split()[0] for gene in sep.split(target.name))
        target.name = f"{target.chrom}:{target.start}-{target.stop}-{target.name}"
        bait2genes[target.name] = genes
        
    details = {"targets": targets,
               "max_fragment_size": max_fragment_size,
               "retain_offtarget": retain_offtarget}
    stats = {"fragments_per_target": {bait.name: 0 for bait in targets},
             "ontarget": 0,
             "offtarget": 0}
    
    
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


    ontarget = stats.pop("ontarget")
    offtarget = stats.pop("offtarget")
    stats["offtarget"] = float(offtarget) / (ontarget + offtarget)
    
    include = set()
    exclude = set()
    for target in cnv.split():
        if target.startswith("-"):
            exclude.add(target[1:])
        else:
            include.add(target)

    if include:
        baseline = []
        depths = defaultdict(list)
        for bait, depth in stats["fragments_per_target"].items():
            genes = bait2genes[bait]
            matches = genes & include
            for gene in matches:
                depths[gene].append(depth)
            if not matches and not genes & exclude:
                baseline.append(depth)
        
        baseline = mean(baseline)
        stats["copies_per_cell"] = {t: mean(d) / baseline for t, d in depths.items()}

    save_stats(stats_file, stats)



def _filter_read(read, stats, targets, max_fragment_size, retain_offtarget):
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
        sys.exit("SAM file not filtered by name or corrupt")
    
    # Both mapped to same reference and pointing in opposite directions
    # therefore may be a concordant read pair.
    if not (primary[0][FLAG] & UNMAPPED) and not (primary[1][FLAG] & UNMAPPED) and \
        primary[0][RNAME] == primary[1][RNAME] and \
        primary[0][FLAG] & RC != primary[1][FLAG] & RC:
        
        if primary[0][FLAG] & RC:
            primary = primary[::-1]
        
        left_cigar = string2cigar(primary[0][CIGAR])
        right_cigar = string2cigar(primary[1][CIGAR])
        start = int(primary[0][POS])
        stop = int(primary[1][POS]) + cigar_len(right_cigar, CONSUMES_REF) - 1
        size = stop - start + 1
        for num, op in chain(left_cigar, right_cigar):
            if op in "IS":
                size += num
            elif op in "DN":
                size -= num
        if size < 0 or (max_fragment_size and size > max_fragment_size):
            size = 0
        
    else:
        size = 0
    

    if size:
        segments = [Entry(primary[0][RNAME], start, stop)]
    else:
        segments = []
        non_primary = read
    for segment in non_primary:
        start = int(segment[POS])
        segments.append(Entry(segment[RNAME], start, start + cigar_len(string2cigar(segment[CIGAR]), CONSUMES_REF)))
    
    
    match = None
    for segment in segments:
        for target in targets.touched_by(segment):
            offset = abs(segment.start + segment.stop - target.start - target.stop)
            try:
                if offset > best_offset:
                    continue
            except NameError:
                pass
            best_offset = offset
            match = target.name
    
    if match:
        stats["fragments_per_target"][match] += 1
        stats["ontarget"] += 1
    else:
        stats["offtarget"] += 1
        if not retain_offtarget:
            return ""

    for segment in read:
        segment[FLAG] = str(segment[FLAG])
    return "".join("\t".join(segment) for segment in read)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', help="Input sam file, must be sorted by name.")
    parser.add_argument("-b", "--bed", help="Bed file of on-target regions.", dest="bed_file", required=True)
    parser.add_argument("-o", "--output", help="Output file.", dest="output_file", default=argparse.SUPPRESS)
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="stats_file", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--max-fragment-size", help="Maximum fragment size to be considered a genuine read pair.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-r", "--retain-offtarget", help="Retain offtarget reads in output.", action='store_const', const=True)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-c", "--cnv", help="Target names over which to calculate copy numbers. " \
                                                "names preceeded by - will be excluded from baseline.", default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        ontarget(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

