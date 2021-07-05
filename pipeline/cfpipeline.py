#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob

from pipeline import run, Pipe, guess_sample_name



def cfpipeline():
    """Cell free pipeline.
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Paths of input fastq or fastq.gz files. Order is important if paired end reads.")
    parser.add_argument("-r", "--reference", help="Path to reference genome or containing directory.", required=True)
    parser.add_argument("-n", "--name", help="Sample name used to name output files. Will be guessed from input fastq if not provided", default="")
    parser.add_argument("-p", "--panel", help="Path to covermi panel which must contain targets bedfile.", default="")
    parser.add_argument("-u", "--umi", help="UMI type (prism, thruplex_hv or thruplex) or empty strng if no umis.", default="")
    parser.add_argument("-v", "--vep", help="Path to vep data.", default="")
    parser.add_argument("-m", "--min-family-size", help="Minimum family size. Families smaller than this will be filtered", type=int, default=1)
    parser.add_argument("-f", "--min-vaf", help="Minimum variant allele frequency for a variant to be called when using VarDict.", type=float, default=None)
    parser.add_argument("-c", "--cnv", help="Whitespace separated list of target names, as specified in targets bedfile, over which to calculate copy number variation.", default="")
    parser.add_argument("-d", "--sizes", help="Whitespace separated list of reference names over which to calculate fragment size distribution.", default="")
    parser.add_argument("-b", "--translocations", help="Call translocations (supplementary reads aligned to different chromosomes).", action="store_const", const=True, default=False)
    parser.add_argument("-i", "--interleaved", help="Each input fastq contains alternating reads 1 and 2.", action="store_const", const=True, default=False)
    parser.add_argument("-o", "--output", help="Path to write output files to.", default=".")
    parser.add_argument("-t", "--threads", help="Number of threads to use, defaults to all available threads if not specified.", type=int, default=None)
    a = parser.parse_args()
    
    threads = a.threads or run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    
    if not a.name:
        a.name = guess_sample_name(a.input_fastqs)
        if not a.name:
            sys.exit("Ambiguous sample name")
    
    if a.min_vaf is None:
        a.min_vaf = 0.01 if a.min_family_size == 1 else 0.001
    
    a.reference = os.path.abspath(a.reference)
    a.input_fastqs = [os.path.abspath(path) for path in a.input_fastqs]
    if a.panel:
        a.panel = os.path.abspath(a.panel)
    if a.vep:
        a.vep = os.path.abspath(a.vep)
    os.chdir(a.output)
    
    a.reference = (glob.glob(f"{a.reference}/*.fna") + [a.reference])[0]
    targets_bedfile = (glob.glob(f"{a.panel}/*.bed") + [None])[0] if a.panel else ""
    stats = f"{a.name}.stats.json"
    pipe = Pipe()
    
    
    # Remove umis and do some basic fastq qc
    interleaved_fastq = f"{a.name}.interleaved.fastq"
    command = ["udini", "--output", interleaved_fastq,
                        "--stats", stats,
                        "--umi", a.umi]
    if a.interleaved:
        command.append("--interleaved")
    pipe(command + a.input_fastqs)
    
    
    base_sam = f"{a.name}.base.sam"
    with open(base_sam, "wb") as f_out:
        pipe(["bwa", "mem", "-t", threads, 
                            "-p", # interleaved paired end fastq
                            "-C", # Append fastq comment to sam
                            "-Y", # Soft clip non-primary reads
                            a.reference, 
                            interleaved_fastq], stdout=f_out)
    os.unlink(interleaved_fastq)
    
    
    namesorted_sam = f"{a.name}.namesorted.sam"
    pipe(["samtools", "sort", "-n", # sort by name
                              "-o", namesorted_sam,
                              "-@", threads,
                              base_sam])
    os.unlink(base_sam)
    
    
    trimmed_sam = f"{a.name}.trimmed.sam"
    pipe(["trim_sam", "--output", trimmed_sam,
                      "--stats", stats,
                      "--threads", threads,
                      namesorted_sam])
    os.unlink(namesorted_sam)
    
    
    sorted_sam = f"{a.name}.sorted.sam"
    pipe(["samtools", "sort", "-o", sorted_sam,
                              "-@", threads,
                              trimmed_sam])
    os.unlink(trimmed_sam)
    
    deduplicated_sam = f"{a.name}.deduplicated.sam"
    pipe(["elduderino_py", "--output", deduplicated_sam,
                           "--stats", stats,
                           "--min-family-size", a.min_family_size,
                           "--umi", a.umi,
                           "--threads", threads,
                           sorted_sam])
    os.unlink(sorted_sam)
    
    
    namesorted_sam = f"{a.name}.namesorted.sam"
    pipe(["samtools", "sort", "-n", # sort by name
                              "-o", namesorted_sam,
                              "-@", threads,
                              deduplicated_sam])
    os.unlink(deduplicated_sam)
    
    
    pipe(["size", "--stats", stats,
                  "--rnames", a.sizes,
                  "--output", f"{a.name}.sizes.pdf",
                  namesorted_sam])
    
    
    ontarget_sam = f"{a.name}.ontarget.sam"
    pipe(["ontarget", "--output", ontarget_sam,
                      "--bed", targets_bedfile,
                      "--stats", stats,
                      "--cnv", a.cnv,
                      "--threads", threads,
                      namesorted_sam])
    os.unlink(namesorted_sam)


    namesorted_sam = f"{a.name}.namesorted.sam"
    pipe(["samtools", "sort", "-n", # sort by name
                              "-o", namesorted_sam,
                              "-@", threads, 
                              ontarget_sam])
    os.unlink(ontarget_sam)
    

    fixed_sam = f"{a.name}.fixed.sam"
    pipe(["samtools", "fixmate", namesorted_sam, fixed_sam])
    os.unlink(namesorted_sam)
        
        
    if a.translocations:
        pipe(["breakpoint", "--output", f"{a.name}.translocations.tsv",
                            fixed_sam])
    
    
    bam = f"{a.name}.bam"
    pipe(["samtools", "sort", "-o", bam,
                              "-@", threads, 
                              fixed_sam])
    pipe(["samtools", "index", bam])
    os.unlink(fixed_sam)
    
    
    if a.panel:
        pipe(["covermi_stats", "--panel", a.panel,
                               "--output", f"{a.name}.covermi.pdf",
                               "--stats", stats,
                               bam])


    mpileup = f"{a.name}.mpileup"
    pipe(["samtools", "mpileup", "-o", mpileup,
                                 "-f", a.reference,
                                 "-A",
                                 "-B",
                                 "-q", "10",
                                 "-d", "10000000",
                                 bam])

    pvalue_vcf = f"{a.name}.pvalue.vcf"
    with open(pvalue_vcf, "wb") as f_out:
        pipe(["varscan", "mpileup2cns", mpileup,
                                        "--variants",
                                        "--output-vcf", "1",
                                        "--min-coverage", "1",
                                        "--min-var-freq", "0",
                                        "--min-avg-qual", "20",
                                        "--min-reads2", "3",
                                        "--p-value", "0.05",
                                        "--strand-filter", "1"], stdout=f_out)
    os.unlink(mpileup)
    
    vcf = f"{a.name}.varscan.vcf"
    pipe(["vcf_pvalue_2_phred", pvalue_vcf, "--output", vcf])
    os.unlink(pvalue_vcf)
    
    if a.vep:
        pipe(["annotate_panel", "--vep", a.vep,
                                "--output", f"{a.name}.varscan.annotation.tsv",
                                "--threads", threads,
                                "--panel", a.panel,
                                vcf])
    
    
    vardict_table = f"{a.name}.vardict.tsv"
    with open(vardict_table, "wb") as f_out:
        pipe(["vardictjava", "-K", # include Ns in depth calculation
                             "-deldupvar", # variants are only called if start position is inside the region interest
                             "-G", a.reference,
                             "-N", a.name,
                             "-b", bam,
                             "-Q", "10",
                             "-f", a.min_vaf,
                             "-th", threads,
                             "-u", # count mate pair overlap only once
                             "-fisher", # perform work of teststrandbias.R
                             targets_bedfile], stdout=f_out)
    
    vcf = f"{a.name}.vcf"
    with open(vardict_table, "rb") as f_in:
        with open(vcf, "wb") as f_out:
              pipe(["var2vcf_valid.pl", "-A", # output all variants at same position
                                        "-f", a.min_vaf,
                                        "-N", a.name], stdin=f_in, stdout=f_out)
    os.unlink(vardict_table)
    
    if a.vep:
        pipe(["annotate_panel", "--vep", a.vep,
                                "--output", f"{a.name}.annotation.tsv",
                                "--reference", a.reference,
                                "--threads", threads,
                                "--panel", a.panel,
                                vcf])
    
    
    #vaf_plot = f"{a.name}.vaf.pdf"
    pipe(["vcf_stats", vcf, 
                       "--stats", stats])
                       #"--output", vaf_plot])
    
    print(pipe.durations, file=sys.stderr, flush=True)



def main():
    try:
        cfpipeline()
    except OSError as e:
        # File input/output error. This is not an unexpected error therfore
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

