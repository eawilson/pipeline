import pdb
import argparse
import json
import sys
import os
import glob
from collections import Counter
from multiprocessing import Process, Queue
from collections.abc import Mapping
from itertools import chain

from .utils import run

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
READX = READ1 | READ2

CONSUMES_REF = "MDN=X"
CONSUMES_READ = "MIS=X"

L = 0
R = 1

RCOMPLEMENT = str.maketrans("ATGC", "TACG")



def ltrim(seg, bases):
    pos = int(seg[POS])
    cigar = []
    for num, op in string2cigar(seg[CIGAR]):
        if bases:
            if op in CONSUMES_READ:
                if num < bases:
                    bases -= num
                    if op in CONSUMES_REF:
                        pos += num
                else:
                    if op in CONSUMES_REF:
                        pos += bases
                    num -= bases
                    bases = 0
            
            elif op in CONSUMES_REF:
                pos += num

        if num and not bases:
            cigar.append((num, op))
    seg[CIGAR] = cigar2string(cigar)
    seg[POS] = str(pos)



def rtrim(seg, bases):
    cigar = string2cigar(seg[CIGAR])
    while bases and cigar:
        try:
            num, op = cigar.pop()
        except IndexError:
            sys.exit("Malformed cigar string {}".format())
        if op in CONSUMES_READ:
            if num <= bases:
                bases -= num
            else:
                cigar.append((num - bases, op))
                bases = 0
    seg[CIGAR] = cigar2string(cigar)



def worker_process(input_queue, output_queue, stats_queue, stats, details):
    while True:
        read = input_queue.get()
        if read is None:
            break
        output = _trim_read(read, stats, details)
        if output:
            output_queue.put(output)
    


def filewriter_process(output_queue, output_file):
    with open(output_file, "wt") as f_out:
        while True:
            output = output_queue.get()
            if output is None:
                break
            f_out.write(output)



def string2cigar(cigstr):
    if cigstr == "*":
        return []
    
    cig = []
    num = ""
    for char in cigstr:
        if char.isnumeric():
            num += char
        else:
            try:
                cig.append((int(num), char))
            except ValueError:
                sys.exit(f"Malformed cigar string {cigstr}")
            num = ""
    if num:
        raise sys.exit(f"Malformed cigar string {cigstr}")
    return cig



def cigar2string(cig):
    return "".join(str(val) for val in chain(*cig)) or "*"



def trim_sam(input_sam, output_file="output.trimmed.sam", threads=0, stats_file="stats.json"):

    threads = 1
    if not threads:
        threads = int(run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip())
    
    details = {}
    stats = {"overlap": 0,
             "mismatches": 0}

    multithreaded = (threads > 1)
    with (open(output_file, "wt") if not multithreaded else nullcontext()) as f_out:
        
        if multithreaded:
            input_queue = Queue()
            output_queue = Queue()
            stats_queue = Queue()
            workers = []
            for i in range(threads - 1):
                workers.append(Process(target=worker_process, args=(input_queue, output_queue, stats_queue, stats, details)))
            workers.append(Process(target=filewriter_process, args=(output_queue, output_file)))
            for worker in workers:
                worker.start()
            trim_read = input_queue.put
            write = output_queue.put
            
        else:
            def trim_read(read):
                output = _trim_read(read, stats, details)
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
                
                qname = row[:row.index("\t")]
                if qname != current_qname:
                    if multithreaded and not all(worker.is_alive() for worker in workers):
                        sys.exit("Worker thread unexpectedly terminated")
                    trim_read(read)
                    current_qname = qname
                    read = []
                read.append(row)
                
            trim_read(read)


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


    try:
        with open(stats_file, "rt") as f:
            old_stats = json.load(f)
    except OSError:
        old_stats = {}
    old_stats["sequencing_error_rate"] = float(stats["mismatches"]) / stats["overlap"]
    with open(stats_file, "wt") as f:
        json.dump(old_stats, f, sort_keys=True, indent=4)



def _trim_read(read, stats, details):
    if not read:
        return ""
    
    primary = []
    segments = [segment.split("\t") for segment in read]
    for segment in segments:
        segment[FLAG] = int(segment[FLAG])
        if not segment[FLAG] & SEC_OR_SUP: # Primary read
            primary.append(segment)

    if len(primary) != 2:
        sys.exit("Multiple primary reads in SAM file")
    
    # Unmapped, mappend to a different reference or pointing in the same direction
    # therefore cannot be a concordant pair
    if primary[L][FLAG] & UNMAPPED or primary[R][FLAG] & UNMAPPED or primary[L][RNAME] != primary[R][RNAME] or primary[L][FLAG] & RC == primary[R][FLAG] & RC:
        return "".join(read)
        
    
    # Ensure the left segment has the lowest ref pos
    lref = int(primary[L][POS])
    rref = int(primary[R][POS])
    if rref < lref:
        primary = primary[::-1]
        lref, rref = rref, lref
        
    lref -= 1
    lread = -1
    for num, op in string2cigar(primary[L][CIGAR]):
        if op in CONSUMES_REF:
            if lref + num > rref:
                num = rref - lref
            lref += num
        if op in CONSUMES_READ:
            lread += num
        
        if lref == rref:
            for num, op in string2cigar(primary[R][CIGAR]):
                if op in CONSUMES_REF:
                    break
                if op in CONSUMES_READ:
                    lread -= num
            break
    
    # Segments don't touch therefore return
    else:
        return "".join(read)
    
    
    # Now ensure the left segment is correctly orientated
    if primary[L][FLAG] & RC:
        primary = primary[::-1]
        lread = -lread
    
    # lread is now the position of the base in the left read that
    # overlaps the first base in the right read. If lread is less
    # than zero then there is readthrough into the opposite umi
    
    lseq = list(primary[L][SEQ])
    lqual = list(primary[L][QUAL])
    rseq = list(primary[R][SEQ])
    rqual = list(primary[R][QUAL])
    
    # Overhang at beginning of right read
    roverhang = max(-lread, 0)
    if roverhang:        
        rseq = rseq[roverhang:]
        rqual = rqual[roverhang:]
        
    # Overhang at end of left read
    loverhang = max(len(primary[L][SEQ]) - lread - len(primary[R][SEQ]), 0)
    if loverhang:
        lseq = lseq[:-loverhang]
        # QUAL will end in \n if no optional fields afterwards
        if lqual[-1] == "\n":
            lqual[-loverhang-1] = "\n"
        lqual = lqual[:-loverhang]
    
    mismatches = 0
    offset = max(lread, 0)
    overlap = len(lseq) - offset
    for i in range(0, overlap):
        if lseq[i+offset] != rseq[i]:
            mismatches += 1
            if ord(lqual[i+offset]) > ord(rqual[i]) + 10:
                rseq[i] = lseq[i+offset]
                rqual[i] = lqual[i+offset]
            elif ord(rqual[i]) > ord(lqual[i+offset]) + 10:
                lseq[i+offset] = rseq[i]
                lqual[i+offset] = rqual[i]
            else:
                rseq[i] = "N"
                rqual[i] = "!"
                lseq[i+offset] = "N"
                lqual[i+offset] = "!"

    seq = (["".join(lseq)], ["".join(rseq)])
    qual = (["".join(lqual)], ["".join(rqual)])
    lbases = ((0, loverhang), (roverhang, 0))
    rbases = ((loverhang, 0), (0, roverhang))
    
    stats["overlap"] += overlap
    stats["mismatches"] += mismatches
    
    for segment in segments:
        # r1r2 = L if corresponds to left primary segment and R if corresponds to right primary segment
        r1r2 = segment[FLAG] & READX == primary[R][FLAG] & READX
        # rc = 0 if orientated] the same way as the corresponding primaary segment and 1 if reversed
        rc = segment[FLAG] & RC != primary[r1r2][FLAG] & RC
        
        bases = lbases[r1r2][rc]
        if bases:
            ltrim(segment, bases)
            
        else:
            bases = rbases[r1r2][rc]
            if bases:
                rtrim(segment, bases)
        
        try:
            segment[SEQ] = seq[r1r2][rc]
        except IndexError:
            seq[r1r2].append(seq[r1r2][0][::-1].translate(RCOMPLEMENT))
            segment[SEQ] = seq[r1r2][rc]
            if qual[r1r2][0][-1] != "\n":
                qual[r1r2].append(qual[r1r2][0][::-1])
            else:
                qual[r1r2].append("{}\n".format(qual[r1r2][0][len(qual[r1r2][0])-2::-1]))
        segment[QUAL] = qual[r1r2][rc]
        
    for segment in segments:
        segment[FLAG] = str(segment[FLAG])
    
    return "".join("\t".join(segment) for segment in segments)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', help="Input sam file, must be sorted by name.")
    parser.add_argument("-o", "--output", help="Output sam file.", dest="output_file", default=argparse.SUPPRESS)
    parser.add_argument("-s", "--stats", help="Statistics file.", dest="stats_file", default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        trim_sam(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

