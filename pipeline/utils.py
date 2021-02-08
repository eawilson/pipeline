import subprocess
import os
import sys
import re
import json
import pdb
import datetime
import shlex
from collections import defaultdict, Counter
from collections.abc import Mapping
from itertools import chain



__all__ = ["run", "pipe", "Pipe", "save_stats", "string2cigar", "cigar2string", "guess_sample_name", "nullcontext", "CONSUMES_REF", "CONSUMES_READ"]


CONSUMES_REF = "MDN=X"
CONSUMES_READ = "MIS=X"



class nullcontext(object):
    def __enter__(self):
        return None

    def __exit__(self, *excinfo):
        pass

#ILLUMINA_FASTQ = re.compile(r"(.+)_S([0-9]{1,2})_L([0-9]{3})_R([12])_001\.fastq(\.gz)?$") # name, s_number, lane, read, gzip



#def illumina_readgroup(filepath):
    #basename = os.path.basename(filepath)
    #sample = "_".join(basename.split("_")[:-4]) # remove the _Sx_Lxxx_Rx_001.fastq from the name
    #with open(filepath) as f:
        #identifier = f.readline().split(":")
    #flowcell = identifier[2]
    #return "@RG\\tID:{}\\tSM:{}".format(flowcell, sample)



def guess_sample_name(fastqs):
    ILLUMINA_FASTQ = re.compile(r"(_S[0-9]{1,2}_L[0-9]{3}_R[12]_001)?\.fastq(\.gz)?$")
    names = set()
    for fastq in fastqs:
        match = ILLUMINA_FASTQ.search(fastq)
        if match:
            fastq = fastq[:match.start()]
        names.add(fastq)
    if len(names) == 1:
        return list(names)[0]    



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



def rekey(mapping):
    """ Recursively convert all numeric text keys to integer keys. This will
        enable correct ordering when re-written to file.
    """
    new = {}
    for k, v in mapping.items():
        try:
            k = int(k)
        except ValueError:
            pass
        if isinstance(v, Mapping):
            v = rekey(v)
        new[k] = v
    return new



def save_stats(path, update):
    try:
        with open(path, "rt") as f_in:
            stats = rekey(json.load(f_in))
        stats.update(update)
    except OSError:
        stats = update
    
    with open(path, "wt") as f_out:
        json.dump(stats, f_out, sort_keys=True, indent=4)



def pretty_duration(seconds):
    mins, secs = divmod(int(seconds), 60)
    hours, mins = divmod(mins, 60)
    duration = [f"{hours} hours"] if hours else []
    if hours or mins:
        duration.append(f"{mins} minutes")
    duration.append(f"{secs} seconds")
    return " ".join(duration) 



class Pipe(object):
    """ Wrapper arond the pipe function that will maintain a record of the
        time taken to run each command. This is stored by command, ie if
        a single command is run several times the time will be recorded as
        the total time of all of the invocations.
    """
    def __init__(self):
        self._durations = Counter()
        
    def __call__(self, *args, **kwargs):
        start = datetime.datetime.now()
        ret = pipe(*args, **kwargs)
        stop = datetime.datetime.now()
        self._durations[args[0][0]] += (stop - start).total_seconds()
        return ret
    
    @property
    def durations(self):
        padding = max(len(key) for key in self._durations)
        template = f"{{:{padding}}} {{}} "
        return "\n".join(template.format(k, pretty_duration(v)) for k, v in sorted(self._durations.items()))



def pipe(args, exit_on_failure=True, **kwargs):
    """ Runs a main pipeline command. Output is bytes rather than string and
        is expected to be captured via stdout redirection or ignored if not
        needed. The command is echoed to stderr before the command is run.
    """
    args = [str(arg) for arg in args]
    print(" ".join(shlex.quote(arg) for arg in args), file=sys.stderr, flush=True)
    completedprocess = subprocess.run(args, **kwargs)
    sys.stderr.flush()
    if exit_on_failure and completedprocess.returncode:
        sys.exit(completedprocess.returncode)
    return completedprocess



def run(args, exit_on_failure=True):
    """ Run a unix command as a subprocess. Stdout and stderr are captured as
        a string for review if needed. Not to be used for main pipeline
        comands which should be called with pipe instead.  
    """
    args = [str(arg) for arg in args]
    completedprocess = subprocess.run(args,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      universal_newlines=True)
    if exit_on_failure and completedprocess.returncode:
        for line in completedprocess.stderr.splitlines():
            print(line, file=sys.stderr, flush=True)
        sys.exit(completedprocess.returncode)
    return completedprocess
