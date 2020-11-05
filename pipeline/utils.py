import subprocess
import os
import sys
import re
import pdb
import datetime
from collections import defaultdict, Counter

import covermi
from boto3 import client



__all__ = ["run", "pipe", "Pipe"]

#ILLUMINA_FASTQ = re.compile(r"(.+)_S([0-9]{1,2})_L([0-9]{3})_R([12])_001\.fastq(\.gz)?$") # name, s_number, lane, read, gzip



#def illumina_readgroup(filepath):
    #basename = os.path.basename(filepath)
    #sample = "_".join(basename.split("_")[:-4]) # remove the _Sx_Lxxx_Rx_001.fastq from the name
    #with open(filepath) as f:
        #identifier = f.readline().split(":")
    #flowcell = identifier[2]
    #return "@RG\\tID:{}\\tSM:{}".format(flowcell, sample)



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
    quoted_args = [(f'"{arg}"' if " " in arg else arg) for arg in args]
    print(" ".join(quoted_args), file=sys.stderr, flush=True)
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


