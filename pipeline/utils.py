import subprocess
import os
import sys
import pdb
import re
import csv
import io
import argparse
from itertools import zip_longest
import pdb
from collections import defaultdict, namedtuple

import covermi
from boto3 import client



__all__ = ["run", "mount_basespace", "unmount_basespace", "mount_instance_storage", "list_basespace_fastqs", "ungzip_and_combine_illumina_fastqs", \
            "load_panel_from_s3", "illumina_readgroup", "pipe", "command_line_arguments", "sample_name"]

BUCKET = "omdc-data"
ILLUMINA_FASTQ = re.compile(r"(.+)_S([0-9]{1,2})_L([0-9]{3})_R([12])_001\.fastq(\.gz)?$") # name, s_number, lane, read, gzip


def command_line_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastqs', nargs="+", help="Fastq files.")
    parser.add_argument("-p", "--panel", help="Directory containing panel data.")
    parser.add_argument("-g", "--genome", help="Reference genome.")
    parser.add_argument("-b", "--bam", help="Matched normal bam to perform tumour/normal subtraction.", default=argparse.SUPPRESS)
    parser.add_argument("-o", "--output", help="Output directory.", dest="dest_dir", default=argparse.SUPPRESS)
    args = parser.parse_args()
    return vars(args)



def sample_name(fastqs):
    """ Check all fastqs have an illumina format file name and all belong to
        the same sample. Return the sample name.
        Args:
            fastqs (list of str): List of fastq paths
            
        Returns:
            str: Sample name.
            
        Raises:
            RuntimeError: If sample name is not consistent between samples.
    """
    matches = [ILLUMINA_FASTQ.match(os.path.basename(fastq)) for fastq in fastqs]
    if not None in matches:
        names = set(match.group(1) for match in matches)
        if len(names) == 1:
            return list(names)[0]
    raise RuntimeError("Inconsistent fastq names.")
    
    

def fasta_path(path):
    """ Return path to reference genome fasta which may be located at an
        arbitrary depth depth within path.
        
        Args:
            path (str): Path to fasta file or a containing directory.
            
        Returns:
            str: Path to reference fasta file.
            
        Raises:
            RuntimeError if no fasta is found or if multiple fastas are found.
    """
    if os.path.isfile(path):
        return path
    fastas = []
    for root, dns, fns in os.walk(path):
        for fn in fns:
            if fn.endswith(".fna"):
                fastas += [os.path.join(root, fn)]
    if len(fastas) != 1:
        raise RuntimeError(f"{path} must contain a single fasta file.")
    return fastas[0]



def ungzip_and_combine_illumina_fastqs(input_fastqs, sample, overwrite=False):
    """ Accepts list of fastqs to be ungzipped (if needed) and combined.
        Output will be written to current working directory.
        
        Args:
            inoput_fastqs (list of str): List of fastq paths.
            sample (str): Sample name.
            overwrite (bool): If False then don't overwrite fastqs if exist.
            
        Returns:
            list of str: List of ungzipped merged fastq paths or None if the
                input_fastqs already consist of two ungzipped fastqs.
            
        Raises:
            RuntimeError if input_fastqs are not paired.
    """
    input_fastqs = [os.path.abspath(fastq) for fastq in input_fastqs]
    input_fastqs = sorted(input_fastqs, key=lambda p:os.path.basename(p))
    r1_fastqs = input_fastqs[0::2]
    r2_fastqs = input_fastqs[1::2]
    
    # Check that fastqs are paired before we merge them.
    for r1_fastq, r2_fastq in zip_longest(r1_fastqs, r2_fastqs, fillvalue=""):
        for r1_char, r2_char in zip_longest(os.path.basename(r1_fastq),
                                            os.path.basename(r1_fastq),
                                            fillvalue=""):
            if r1_char != r2_char and (r1_char != "1" or r2_char != "2"):
                raise RuntimeError(f"{r1_fastq} and {r2_fastq} are not paired.")
            
    # If we only have a single pair of ungzipped fastqs then we are done.
    if len(r1_fastqs) == 1 and r1_fastqs[0].endswith(".fastq"):
        return None
    
    ret = []
    for read, fastqs in (("1", r1_fastqs), ("2", r2_fastqs)):
        output_fastq = f"{sample}_R{read}.fastq"
        ret += [output_fastq]
        if overwrite or not os.path.exists(output_fastq):
            with open(output_fastq, "wb") as f_out:
                for fastq in fastqs:
                    if fastq.endswith(".gz"):
                        subprocess.run(["gzip", "-dc", fastq], stdout=f_out)      
                    else:
                        subprocess.run(["cat", fastq], stdout=f_out)
    return ret



def illumina_readgroup(filepath):
    basename = os.path.basename(filepath)
    sample = "_".join(basename.split("_")[:-4]) # remove the _Sx_Lxxx_Rx_001.fastq from the name
    with open(filepath) as f:
        identifier = f.readline().split(":")
    flowcell = identifier[2]
    return "@RG\\tID:{}\\tSM:{}".format(flowcell, sample)



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



def pipe(args, exit_on_failure=True, **kwargs):
    """ Runs a main pipeline command. Output is bytes rather than string and
        is expected to be captured via stdout redirection or ignored if not
        needed. The command is echoed to stderr before the command is run.
    """
    args = [str(arg) for arg in args]
    print(" ".join(args), file=sys.stderr, flush=True)
    completedprocess = subprocess.run(args, **kwargs)
    sys.stderr.flush()
    if exit_on_failure and completedprocess.returncode:
        sys.exit(completedprocess.returncode)



def mount_basespace():
    basespace_dir = os.path.join(os.path.expanduser("~"), "basespace")
    if not os.path.exists(basespace_dir):
        os.mkdir(basespace_dir)
        
    completed = run(["mount"])
    for line in completed.stdout.split("\n"):
        if line.startswith("basemount on {}".format(basespace_dir)):
            print("Basespace already mounted.")
            return basespace_dir
        
    print("Mounting basespace.")
    run(["basemount", basespace_dir])
    return basespace_dir



def unmount_basespace():
    basespace_dir = os.path.join(os.path.expanduser("~"), "basespace")
    completed = run(["mount"])
    for line in completed.stdout.split("\n"):
        if line.startswith("basemount on {}".format(basespace_dir)):
            print("Unmounting basespace.")
            run(["basemount", "--unmount", basespace_dir])
            return
        
    print("Basespace not mounted.")



def list_basespace_fastqs(project="", sample=""):
    basespace_path = os.path.join(os.path.expanduser("~"), "basespace")
    project_regex = re.compile(project)
    sample_regex = re.compile(sample)

    matches = []
    projects_dir = os.path.join(basespace_path, "Projects")
    projects = [project for project in os.listdir(projects_dir) if not project.startswith(".") and project_regex.search(project)]    
    
    for project in projects:
        print("Indexing basespace project {}.".format(project))
        samples_dir = os.path.join(projects_dir, project, "Samples")
        samples = [sample for sample in os.listdir(samples_dir) if not sample.startswith(".") and sample_regex.search(sample)]
        
        for sample in samples:
            files_dir = os.path.join(samples_dir, sample, "Files")
            matches += [os.path.join(files_dir, fn) for fn in os.listdir(files_dir) if not fn.startswith(".") and (fn.endswith(".fastq") or fn.endswith(".fastq.gz"))]
            
    return matches



