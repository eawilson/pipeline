#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob

from pipeline import run, pipe, vcf_pvalue_2_phred, create_report, Pipe
from covermi import Panel, covermimain



def cfpipeline(sample, input_fastqs, reference, panel, umi=None, vep=None, min_family_size=1, max_fragment_size=None, threads=None):
    """Cell free pipeline.

    Args:
        sample (str): Sample name, used to name output files and covermi plot.
        input_fastqs (list of str): Paths of input paired fastqs or fastq.gzs,
            If paired fastqs then order of files is important.
        reference (str): Path to reference fasta or containing directory.
        panel (str): Path to covermi panel which must contain targets bedfile.
        umi (str): umi type or None if no umis.
        vep (str): Path to vep data.
        min_family_size (int): Minimum family size, families smaller
            than this will be filtered.
        max_fragment_size (int): Maximum template length of aligned pair that
            will be considered genuine, pairs with a tlen greater than this
            will be treated as discordant..
        threads (int): Number of threads to use, defaults to all available
            threads if not specified.
        
    Returns:
        None.
    """
    
    if threads is None:
        threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    
    reference = (glob.glob(f"{reference}/*.fna") + [reference])[0]
    targets_bedfile = (glob.glob(f"{panel}/*.bed") + [panel])[0]
    stats = f"{sample}.stats.json"
    pipe = Pipe()
    
    
    # Remove umis and do some basic fastq qc
    interleaved_fastq = f"{sample}.interleaved.fastq"
    udini_options = ["--output", interleaved_fastq,
                     "--stats", stats]
    if umi is not None:
        udini_options += ["--umi", umi]
    pipe(["udini"] + input_fastqs + udini_options)
    
    
    undeduped_unsorted_sam = f"{sample}.undeduped.unsorted.sam"
    with open(undeduped_unsorted_sam, "wb") as f_out:
        pipe(["bwa", "mem", "-t", threads, 
                            "-p", 
                            "-C", 
                            reference, 
                            interleaved_fastq], stdout=f_out)
    os.unlink(interleaved_fastq)
    
    
    undeduped_sam = f"{sample}.undeduped.sam"
    pipe(["samtools", "sort", "-o", undeduped_sam,
                              "-@", threads,
                              undeduped_unsorted_sam])
    os.unlink(undeduped_unsorted_sam)
    
    
    unsorted_sam = f"{sample}.unsorted.sam"
    elduderino_options = ["--output", unsorted_sam,
                          "--stats", stats,
                          "--min-family-size", min_family_size,
                          "--threads", threads]
    if max_fragment_size is not None:
        elduderino_options += ["--max-fragment-size", max_fragment_size]
    if umi is not None:
        elduderino_options += ["--umi", umi]
    if targets_bedfile is not None:
        elduderino_options += ["--bed", targets_bedfile]
    pipe(["elduderino2", undeduped_sam] + elduderino_options)
    os.unlink(undeduped_sam)


    namesorted_unfixed_sam = f"{sample}.namesorted.unfixed.sam"
    pipe(["samtools", "sort", "-n",
                              "-o", namesorted_unfixed_sam,
                              "-@", threads, 
                              unsorted_sam])
    os.unlink(unsorted_sam)
    

    namesorted_sam = f"{sample}.namesorted.sam"
    pipe(["samtools", "fixmate", namesorted_unfixed_sam, namesorted_sam])
    os.unlink(namesorted_unfixed_sam)


    sam = f"{sample}.sam"
    pipe(["samtools", "sort", "-o", sam,
                              "-@", threads, 
                              namesorted_sam])
    os.unlink(namesorted_sam)


    bam = f"{sample}.bam"
    pipe(["samtools", "view", "-b", sam, "-o", bam])
    os.unlink(sam)
    pipe(["samtools", "index", bam])
    

    mpileup = f"{sample}.mpileup"
    mpileup_options = ["-A",
                       "-B",
                       "-q", "10",
                       "-d", "10000000"]
    pipe(["samtools", "mpileup", "-o", mpileup,
                                 "-f", reference] + 
                                 mpileup_options + [bam])
    pvalue_vcf = f"{sample}.pvalue.vcf"
    varscan_options = ["--min-coverage", "1",
                       "--min-var-freq", "0", 
                       "--min-avg-qual", "20",
                       "--min-reads2", "3",
                       "--p-value", "0.05",
                    #"--min-coverage-normal", "1",
                    #"--min-coverage-tumor", "1",
                    #"--min-freq-for-hom","0.75",
                    #"--somatic-p-value", "0.05",
                       "--strand-filter", "1",]
                    #"--validation", "1"


    with open(pvalue_vcf, "wb") as f_out:
        pipe(["varscan", "mpileup2cns", mpileup, "--variants", 
                                                 "--output-vcf", "1"] + 
                                                 varscan_options, stdout=f_out)
    os.unlink(mpileup)
    vcf = f"{sample}.vcf"
    vcf_pvalue_2_phred(pvalue_vcf, vcf)
    os.unlink(pvalue_vcf)
    
    if vep is not None:
        vepjson = f"{sample}.vep"
        vep_options = ["--no_stats",
                    "--dir", vep,
                    "--format", "vcf",
                    "--fork", threads,
                    "--json",
                    "--offline",
                    "--everything",
                    "--force_overwrite"]
        if "refseq" in vep:
            vep_options += ["--refseq"]
        
        pipe(["vep", "-i", vcf,
                     "-o", vepjson,
                     "--fasta", reference] + vep_options)
        create_report(vepjson, panel)
        os.unlink(vepjson)
    
    pipe(["covermi_stats", bam, "--panel", panel,
                                "--output", f"{sample}.covermi.pdf",
                                "--stats", stats,
                                "--sample", sample])
    
    print(pipe.durations, file=sys.stderr, flush=True)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Fastq files.")
    parser.add_argument("-n", "--sample", help="Sample name.", required=True)
    parser.add_argument("-r", "--reference", help="Reference genome.", required=True)
    parser.add_argument("-p", "--panel", help="Directory containing panel data.", required=True)
    parser.add_argument("-u", "--umi", help="Umi.", default=argparse.SUPPRESS)
    parser.add_argument("-v", "--vep", help="Directory containing vep data.")
    parser.add_argument("-m", "--min-family-size", help="Minimum family size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-f", "--max-fragment-size", help="Maximum template legth to be considered a genuine read pair.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    cfpipeline(**vars(args))



if __name__ == "__main__":
    main()

