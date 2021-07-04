#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob
import json


from pipeline import run, Pipe
import subprocess


def multiplexing():
    """Cell free pipeline.
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input_sam', help="Paths of input sorted undeduplicated sam file.")
    parser.add_argument("-n", "--name", help="Sample name used to name output files. Will be guessed from input sam if not provide.", default="")
    parser.add_argument("-u", "--umi", help="UMI type (prism, thruplex_hv or thruplex) or empty strng if no umis.", default="")
    parser.add_argument("-m", "--min-family-size", help="Minimum family size. Families smaller than this will be filtered", type=int, default=1)
    parser.add_argument("-l", "--interval", help="Step size to increment downsampling by.", type=int, required=True)
    parser.add_argument("-r", "--reference", help="Path to reference genome or containing directory.", required=True)
    parser.add_argument("-p", "--panel", help="Path to covermi panel which must contain targets bedfile.", required=True)
    parser.add_argument("-o", "--output", help="Path to write output files to.", default=".")
    args = parser.parse_args()
    
    threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    
    if not args.name:
        args.name = input_sam .split("/")[-1].split(".")[0]
    
    args.reference = os.path.abspath(args.reference)
    args.input_sam = os.path.abspath(args.input_sam)
    args.panel = os.path.abspath(args.panel)
    os.chdir(args.output)
    
    args.reference = (glob.glob(f"{args.reference}/*.fna") + glob.glob(f"{args.reference}/*.fa") + glob.glob(f"{args.reference}/*.fasta") + [args.reference])[0]
    ref_dir = os.path.dirname(args.reference)
    if glob.glob(f"{ref_dir}/*.sa"):
        bwa = "bwa"
    elif glob.glob(f"{ref_dir}/*.0123"):
        bwa = "bwa-mem2"
    else:
        sys.exit("Invalid bwa indexes")
    targets_bedfile = (glob.glob(f"{args.panel}/*.bed") + [None])[0]
    stats = f"{args.name}.stats.json"
    pipe = Pipe()
    
    
    cp = subprocess.run(f"wc -l '{args.input_sam}'", shell=True, stdout=subprocess.PIPE)
    total_reads = int(cp.stdout.encode())
    
    output_file = f"{name}.multiplexing.tsv"
    with open(output_file, "wt") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["sample", "reads", "mean_depth", "mean_family_size", "singleton_rate", "triplicate_plus_rate"])
        
        selected_reads = 0
        while True:
            selected_reads += args.interval
            if selected_reads > total_reads:
                break
            
            
            downsampled_sam = f"{args.name}.downsampled.sam"
            pipe(["downsample_sam", "--output", downsampled_sam,
                                "--number", selected_reads,
                                args.input_sam])
            
            
            deduplicated_fastq = f"{args.name}.deduplicated.fastq"
            pipe(["elduderino", "--output", deduplicated_fastq,
                                "--stats", stats,
                                "--min-family-size", args.min_family_size,
                                "--umi", args.umi,
                                downsampled_sam])
            os.unlink(downsampled_sam)
            
            
            deduplicated_sam = f"{args.name}.deduplicated.sam"
            with open(deduplicated_sam, "wb") as f_out:
                pipe([bwa, "mem", "-t", threads, 
                                "-p", # interleaved paired end fastq
                                "-C", # Append fastq comment to sam
                                "-Y", # Soft clip non-primary reads
                                args.reference, 
                                deduplicated_fastq], stdout=f_out)
            os.unlink(deduplicated_fastq)


            bam = f"{args.name}.bam"
            pipe(["samtools", "sort", "-o", bam,
                                    "-@", threads, 
                                    deduplicated_sam])
            pipe(["samtools", "index", bam])
            os.unlink(deduplicated_sam)
            
            
            pipe(["covermi_stats", "--panel", args.panel,
                                "--output", f"{args.name}.covermi.pdf",
                                "--stats", stats,
                                bam])
            os.unlink(bam)
            
            
            with open(stats, "rt") as f:
                data = json.load(s)
            os.unlink(stats)
            writer.writerow([args.name,
                             selected_reads,
                             data["coverage"]["mean_depth"],
                             data["mean_family_size"],
                             data["singleton_rate"],
                             data["triplicate_plus_rate"]])
            f_out.flush()
            
            
    print(pipe.durations, file=sys.stderr, flush=True)



def main():
    try:
        multiplexing()
    except OSError as e:
        # File input/output error. This is not an unexpected error therfore
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

