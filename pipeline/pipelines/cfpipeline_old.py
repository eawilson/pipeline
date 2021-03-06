#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob

from pipeline import run, Pipe



def cfpipeline(sample, input_fastqs, reference, panel="", umi="", vep="", min_family_size=1, max_fragment_size=0, cnv="", threads=None, rtrim=0, no_elduderino=False):
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
        cnv (str): Whitespace separated list of target names, as specified in
            panel bedfile, over which to calculate copy number variation.
        threads (int): Number of threads to use, defaults to all available
            threads if not specified.
        rtrim (int): Number of bases to trim from the end of read.
        no_elduderino (bool): Do not run elduderino.
        
    Returns:
        None.
    """
    
    if threads is None:
        threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    
    reference = (glob.glob(f"{reference}/*.fna") + [reference])[0]
    targets_bedfile = (glob.glob(f"{panel}/*.bed") + [None])[0] if panel else ""
    stats = f"{sample}.stats.json"
    pipe = Pipe()
    
    
    # Remove umis and do some basic fastq qc
    interleaved_fastq = f"{sample}.interleaved.fastq"
    pipe(["udini", "--output", interleaved_fastq,
                   "--stats", stats,
                   "--umi", umi,
                   "--rtrim", rtrim] + \
                   input_fastqs)
    
    
    undeduped_unsorted_sam = f"{sample}.undeduped.unsorted.sam"
    with open(undeduped_unsorted_sam, "wb") as f_out:
        pipe(["bwa", "mem", "-t", threads, 
                            "-p", 
                            "-C", 
                            "-Y", 
                            reference, 
                            interleaved_fastq], stdout=f_out)
    os.unlink(interleaved_fastq)
    
    
    undeduped_sam = f"{sample}.undeduped.sam"
    pipe(["samtools", "sort", "-o", undeduped_sam,
                              "-@", threads,
                              undeduped_unsorted_sam])
    os.unlink(undeduped_unsorted_sam)
    

    if not no_elduderino:
        unsorted_sam = f"{sample}.unsorted.sam"
        pipe(["elduderino", "--output", unsorted_sam,
                            "--stats", stats,
                            "--min-family-size", min_family_size,
                            "--threads", threads,
                            "--max-fragment-size", max_fragment_size,
                            "--umi", umi,
                            "--bed", targets_bedfile,
                            undeduped_sam])
        #os.unlink(undeduped_sam)
    else:
        unsorted_sam = undeduped_sam


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
        pipe(["panel_copy_numbers", stats, "--targets", cnv])
    
    
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
    
    print(pipe.durations, file=sys.stderr, flush=True)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Fastq files.")
    parser.add_argument("-n", "--sample", help="Sample name.", required=True)
    parser.add_argument("-r", "--reference", help="Reference genome.", required=True)
    parser.add_argument("-p", "--panel", help="Directory containing panel data.", default=argparse.SUPPRESS)
    parser.add_argument("-u", "--umi", help="Umi.", default=argparse.SUPPRESS)
    parser.add_argument("-v", "--vep", help="Directory containing vep data.", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--min-family-size", help="Minimum family size.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-f", "--max-fragment-size", help="Maximum template legth to be considered a genuine read pair.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-c", "--cnv", help="Targets over which to calculate copy numbers.", default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-i", "--rtrim", help="Trim bases from the end of the read.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-e", "--no-elduderino", help="Don't run elduderino.", action="store_const", const=True, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        cfpipeline(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

