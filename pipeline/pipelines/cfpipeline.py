#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob

from pipeline import run, pipe, vcf_pvalue_2_phred, create_report, Pipe
from covermi import Panel, covermimain



def cfpipeline(sample, input_fastqs, reference, panel, umi=None, vep=None, min_family_size=1, threads=None):
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
    
    # Calculate number of available threads with which to run bwa, samtools
    # and vep. Will use all available threads on the assumption that this
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
    
    
    #deduped_interleaved_fastq = f"{sample}.deduped.interleaved.fastq"
    #dude_options = ["--output", deduped_interleaved_fastq,
                    #"--stats", stats]
    #if umi is not None:
        #dude_options += ["--umi", umi]
    #pipe(["dude"] + input_fastqs + dude_options)
    #os.unlink(interleaved_fastq)

    unsorted_unfixed_undeduped_sam = f"{sample}.unsorted.unfixed.undeduped.sam"
    with open(unsorted_unfixed_undeduped_sam, "wb") as f_out:
        pipe(["bwa", "mem", "-t", threads, 
                            "-p", 
                            "-C", 
                            reference, 
                            interleaved_fastq], stdout=f_out)
    os.unlink(interleaved_fastq)


    namesorted_unfixed_undeduped_sam = f"{sample}.namesorted.unfixed.undeduped.sam"
    pipe(["samtools", "sort", "-n",
                              "-o", namesorted_unfixed_undeduped_sam,
                              "-@", threads, 
                              unsorted_unfixed_undeduped_sam])
    os.unlink(unsorted_unfixed_undeduped_sam)
    

    namesorted_undeduped_sam = f"{sample}.namesorted.undeduped.sam"
    pipe(["samtools", "fixmate", namesorted_unfixed_undeduped_sam, namesorted_undeduped_sam])
    os.unlink(namesorted_unfixed_undeduped_sam)
    

    undeduped_sam = f"{sample}.undeduped.sam"
    pipe(["samtools", "sort", "-o", undeduped_sam,
                              "-@", threads, 
                              namesorted_undeduped_sam])
    os.unlink(namesorted_undeduped_sam)
    sys.exit()

    unsorted_sam = f"{sample}.unsorted.sam"
    elduderino_options = ["--output", unsorted_sam,
                          "--stats", stats,
                          "--min-family-size", min_family_size,
                          "--threads", threads]
    if umi is not None:
        elduderino_options += ["--umi", umi]
    if targets_bedfile is not None:
        elduderino_options += ["--bed", targets_bedfile]
    pipe(["elduderino", undeduped_sam] + elduderino_options)
    os.unlink(undeduped_sam)


    sam = f"{sample}.sam"
    pipe(["samtools", "sort", "-o", sam,
                              "-@", threads, 
                              unsorted_sam])
    os.unlink(unsorted_sam)


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

        prop = Panel(panel).properties
        if prop.get("transcript_source", "refseq") == "refseq":
            vep_options += ["--refseq"]
        if "assembly" in prop:
            vep_options += ["--assembly", prop["assembly"]]
        pipe(["vep", "-i", vcf,
                     "-o", vepjson,
                     "--fasta", reference] + vep_options)
        create_report(vepjson, panel)
        os.unlink(vepjson)

    covermi_dir = covermimain(panel, "", bam_path=bam)
    covermi_file = "{}.covermi.tar.gz".format(sample)
    run(["tar", "cvzf", covermi_file, covermi_dir])
    run(["rm", "-r", covermi_dir])
    
    print(pipe.durations, file=sys.stderr, flush=True)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Fastq files.")
    parser.add_argument("-s", "--sample", help="Sample name.", required=True)
    parser.add_argument("-r", "--reference", help="Reference genome.", required=True)
    parser.add_argument("-p", "--panel", help="Directory containing panel data.", required=True)
    parser.add_argument("-v", "--vep", help="Directory containing vep data.")
    parser.add_argument("-u", "--umi", help="Umi.", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--min-family-size", help="Minimum family size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    cfpipeline(**vars(args))



if __name__ == "__main__":
    main()

