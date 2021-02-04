#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob

from pipeline import run, Pipe, guess_sample_name



def cfpipeline(input_fastqs, reference, sample="", panel="", umi="", vep="", min_family_size=1, cnv="", threads=None, interleaved=False, fragment_size_targets=""):
    """Cell free pipeline.

    Args:
        input_fastqs (list of str): Paths of input fastqs or fastq.gzs,
            If paired fastqs then order of files is important.
        reference (str): Path to reference fasta or containing directory.
        sample (str): Sample name, used to name output files and covermi plot.
        panel (str): Path to covermi panel which must contain targets bedfile.
        umi (str): umi type or empty strng if no umis.
        vep (str): Path to vep data.
        min_family_size (int): Minimum family size, families smaller
            than this will be filtered.
        cnv (str): Whitespace separated list of target names, as specified in
            panel bedfile, over which to calculate copy number variation.
        threads (int): Number of threads to use, defaults to all available
            threads if not specified.
        interleaved (bool): Each input fastq is interleaved, containing
            alternating readss 1 and 2.
        
    Returns:
        None.
    """
    
    if threads is None:
        threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    
    if not sample:
        sample = guess_sample_name(input_fastqs)
        if not sample:
            sys.exit("Ambiguous sample name")
    
    reference = (glob.glob(f"{reference}/*.fna") + [reference])[0]
    targets_bedfile = (glob.glob(f"{panel}/*.bed") + [None])[0] if panel else ""
    stats = f"{sample}.stats.json"
    pipe = Pipe()
    
    
    # Remove umis and do some basic fastq qc
    interleaved_fastq = f"{sample}.interleaved.fastq"
    command = ["udini", "--output", interleaved_fastq,
                        "--stats", stats,
                        "--umi", umi]
    if interleaved:
        command.append("--interleaved")
    pipe(command + input_fastqs)
    
    
    unsorted_sam = f"{sample}.unsorted.sam"
    with open(unsorted_sam, "wb") as f_out:
        pipe(["bwa", "mem", "-t", threads, 
                            "-p", # interleaved paired end fastq
                            "-C", # Append fastq comment to sam
                            "-Y", # Soft clip non-primary reads
                            reference, 
                            interleaved_fastq], stdout=f_out)
    os.unlink(interleaved_fastq)
    
    
    untrimmed_sam = f"{sample}.untrimmed.sam"
    pipe(["samtools", "sort", "-n", # sort by name
                              "-o", untrimmed_sam,
                              "-@", threads,
                              unsorted_sam])
    os.unlink(unsorted_sam)
    
    
    trimmed_sam = f"{sample}.trimmed.sam"
    pipe(["trim_sam", "--output", trimmed_sam,
                      "--stats", stats,
                      "--threads", threads,
                      untrimmed_sam])
    os.unlink(untrimmed_sam)
    
    
    undeduplicated_sam = f"{sample}.undeduplicated.sam"
    pipe(["samtools", "sort", "-o", undeduplicated_sam,
                              "-@", threads,
                              trimmed_sam])
    os.unlink(trimmed_sam)
    
    
    deduplicated_sam = f"{sample}.deduplicated.sam"
    pipe(["elduderino", "--output", deduplicated_sam,
                        "--stats", stats,
                        "--min-family-size", min_family_size,
                        "--umi", umi,
                        "--threads", threads,
                        undeduplicated_sam])
    os.unlink(undeduplicated_sam)
    
    
    unfiltered_sam = f"{sample}.unfiltered.sam"
    pipe(["samtools", "sort", "-n", # sort by name
                              "-o", unfiltered_sam,
                              "-@", threads,
                              deduplicated_sam])
    os.unlink(deduplicated_sam)
    
    
    filtered_sam = f"{sample}.filtered.sam"
    pipe(["filter_sam", "--output", filtered_sam,
                        "--bed", targets_bedfile,
                        "--stats", stats,
                        "--threads", threads,
                        unfiltered_sam])
    os.unlink(unfiltered_sam)


    unfixed_sam = f"{sample}.unfixed.sam"
    pipe(["samtools", "sort", "-n", # sort by name
                              "-o", unfixed_sam,
                              "-@", threads, 
                              filtered_sam])
    os.unlink(filtered_sam)
    

    fixed_sam = f"{sample}.fixed.sam"
    pipe(["samtools", "fixmate", unfixed_sam, fixed_sam])
    os.unlink(unfixed_sam)


    bam = f"{sample}.bam"
    pipe(["samtools", "sort", "-o", bam,
                              "-@", threads, 
                              fixed_sam])
    pipe(["samtools", "index", bam])
    os.unlink(fixed_sam)
    

    mpileup = f"{sample}.mpileup"
    pipe(["samtools", "mpileup", "-o", mpileup,
                                 "-f", reference,
                                 "-A",
                                 "-B",
                                 "-q", "10",
                                 "-d", "10000000",
                                 bam])

    pvalue_vcf = f"{sample}.pvalue.vcf"
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
    
    vcf = f"{sample}.vcf"
    pipe(["vcf_pvalue_2_phred", pvalue_vcf, "--output", vcf])
    os.unlink(pvalue_vcf)
    
    
    if cnv:
        pipe(["panel_copy_number", stats, "--targets", cnv])
    
    
    if vep:
        pipe(["annotate_panel", "--vep", vep,
                                "--output", f"{sample}.annotation.tsv",
                                "--threads", threads,
                                "--panel", panel,
                                vcf])
              
    if panel:
        pipe(["covermi_stats", "--panel", panel,
                               "--output", f"{sample}.covermi.pdf",
                               "--stats", stats,
                               "--sample", sample,
                               bam])
    
    fragment_plot = f"{sample}.fragment.sizes.pdf"
    pipe(["fragment_plot", stats, 
                           "--sample", sample,
                           "--output", fragment_plot])
    
    vaf_plot = f"{sample}.vaf.pdf"
    pipe(["vcf_stats", vcf, 
                       "--stats", stats,
                       "--sample", sample,
                       "--output", vaf_plot])
        
    print(pipe.durations, file=sys.stderr, flush=True)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Fastq files.")
    parser.add_argument("-r", "--reference", help="Reference genome.", required=True)
    parser.add_argument("-n", "--sample", help="Sample name.", default=argparse.SUPPRESS)
    parser.add_argument("-p", "--panel", help="Directory containing panel data.", default=argparse.SUPPRESS)
    parser.add_argument("-u", "--umi", help="Umi.", default=argparse.SUPPRESS)
    parser.add_argument("-v", "--vep", help="Directory containing vep data.", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--min-family-size", help="Minimum family size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-c", "--cnv", help="Targets over which to calculate copy numbers.", default=argparse.SUPPRESS)
    parser.add_argument("-d", "--fragment-size-targets", help="Targets over which to calculate fragment size distribution.", default=argparse.SUPPRESS)
    parser.add_argument("-i", "--interleaved", help="Each input fastq contains alternating reads 1 and 2.", action="store_const", const=True, default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        cfpipeline(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

