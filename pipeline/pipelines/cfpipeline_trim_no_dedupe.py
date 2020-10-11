#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob

from pipeline import run, pipe, vcf_pvalue_2_phred, create_report, Pipe
from covermi import Panel, covermimain



def cfpipeline(sample, input_fastqs, reference, panel, umi=None, threads=None):
    """Cell free pipeline.

    Args:
        sample (str): Sample name, used to name output files.
        input_fastqs (list of str): Paths of input paired fastqs or fastq.gzs,
            If paired fastqs then order of files is important.
        reference (str): Path to reference fasta or containing directory.
        panel (str): Path to covermi panel which must contain targets bedfile.
        umi (str): umi type or None if no umis.
        vep (str): Path to vep data.
        min_family_size (int): Minimum family size, families smaller
            than this will be filtered.
        threads (int): Number of threads to use.
        
    Returns:
        None.
    """
    
    # Calculate number of available threads with which to run bwa, samtools,
    # vep, etc. Will use all available threads on the assumption that this
    # is the only program running on the machine.
    if threads is None:
        threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    
    reference = (glob.glob(f"{reference}/*.fna") + [reference])[0]
    targets_bedfile = (glob.glob(f"{panel}/*.bed") + [panel])[0]
    stats = f"{sample}.stats.json"
    
    pipe = Pipe()
    
    interleaved_fastq = f"{sample}.interleaved.fastq"
    udini_options = ["--output", interleaved_fastq,
                     "--stats", stats]
    if umi is not None:
        udini_options += ["--umi", umi]
    pipe(["udini"] + input_fastqs + udini_options)
    
    
    if umi == "prism":
        umi_len = 8
    elif umi == "thruplex_hv":
        umi_len = 7
    else:
        umi_len = 0
    if umi_len:
        clipped_fastq = f"{sample}.clipped.fastq"
        with open(interleaved_fastq, "rt") as f_in:
            with open(clipped_fastq, "wt") as f_out:
                for i, row in enumerate(f_in):
                    if 1 % 2:
                        row = row.rstrip()[:-umi_len]+"\n"
                    f_out.write(row)
        os.unlink(interleaved_fastq)
        os.rename(clipped_fastq, interleaved_fastq)
                    
    
    unsorted_unfixed_sam = f"{sample}.unsorted.unfixed.sam"
    with open(unsorted_unfixed_sam, "wb") as f_out:
        pipe(["bwa", "mem", "-t", threads, 
                            "-p", 
                            "-C", 
                            reference, 
                            interleaved_fastq], stdout=f_out)
    os.unlink(interleaved_fastq)


    namesorted_unfixed_sam = f"{sample}.namesorted.unfixed.sam"
    pipe(["samtools", "sort", "-n",
                              "-o", namesorted_unfixed_sam,
                              "-@", threads, 
                              unsorted_unfixed_sam])
    os.unlink(unsorted_unfixed_sam)
    

    namesorted_sam = f"{sample}.namesorted.sam"
    pipe(["samtools", "fixmate", namesorted_unfixed_sam, namesorted_sam])
    os.unlink(namesorted_unfixed_sam)
    

    sorted_bam = f"{sample}.sorted.bam"
    pipe(["samtools", "sort", "-o", sorted_bam,
                              "-@", threads,
                              namesorted_sam])
    os.unlink(namesorted_sam)

    bam = f"{sample}.bam"
    metrics = f"{sample}.metrics.txt"
    picard_options = []
    if umi is not None:
        picard_options += ["BARCODE_TAG=RX"]
    pipe(["picard", "MarkDuplicates", "I="+sorted_bam, "O="+bam, "M="+metrics] + picard_options)
    os.unlink(sorted_bam)


    pipe(["samtools", "index", bam])


    covermi_dir = covermimain(panel, "", bam_path=bam)
    covermi_file = "{}.covermi.tar.gz".format(sample)
    run(["tar", "cvzf", covermi_file, covermi_dir])
    run(["rm", "-r", covermi_dir])
    
    print(pipe.durations, file=sys.stderr, flush=True)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Fastq files.")
    parser.add_argument("-n", "--sample", help="Sample name.", required=True)
    parser.add_argument("-r", "--reference", help="Reference genome.", required=True)
    parser.add_argument("-p", "--panel", help="Directory containing panel data.", required=True)
    parser.add_argument("-u", "--umi", help="Umi.", default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    cfpipeline(**vars(args))



if __name__ == "__main__":
    main()

