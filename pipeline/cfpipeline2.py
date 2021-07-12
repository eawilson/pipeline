#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob

from pipeline import run, Pipe, guess_sample_name



def cfpipeline2():
    """Cell free pipeline.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Paths of input fastq or fastq.gz files. Order is important if paired end reads.")
    parser.add_argument("-r", "--reference", help="Path to reference genome or containing directory.", required=True)
    parser.add_argument("-n", "--name", help="Sample name used to name output files. Will be guessed from input fastq if not provided", default="")
    parser.add_argument("-p", "--panel", help="Path to covermi panel which must contain targets bedfile.", default="")
    parser.add_argument("-u", "--umi", help="UMI type (prism, thruplex_hv or thruplex) or empty strng if no umis.", default="")
    parser.add_argument("-v", "--vep", help="Path to vep datargs.", default="")
    parser.add_argument("-m", "--min-family-size", help="Minimum family size. Families smaller than this will be filtered", type=int, default=1)
    parser.add_argument("-f", "--min-vaf", help="Minimum variant allele frequency for a variant to be called when using VarDict.", type=float, default=None)
    parser.add_argument("-a", "--min-alt-reads", help="Minimum number of alt reads for a variant to be called.", type=float, default=2)
    parser.add_argument("-c", "--cnv", help="Whitespace separated list of target names, as specified in targets bedfile, over which to calculate copy number variation.", default="")
    parser.add_argument("-d", "--sizes", help="Whitespace separated list of reference names over which to calculate fragment size distribution.", default="")
    parser.add_argument("-b", "--translocations", help="Call translocations (supplementary reads aligned to different chromosomes).", action="store_const", const=True, default=False)
    parser.add_argument("-i", "--interleaved", help="Each input fastq contains alternating reads 1 and 2.", action="store_const", const=True, default=False)
    parser.add_argument("-o", "--output", help="Path to write output files to.", default=".")
    parser.add_argument("-t", "--threads", help="Number of threads to use, defaults to all available threads if not specified.", type=int, default=None)
    parser.add_argument("-s", "--sam-only", help="Quit after producing initial undeduplicated sam.", action="store_const", const=True, default=False)
    parser.add_argument("-C", "--callers", help="Variant callers to use. Valid values are varscan, vardict and mutect2. Defaults to 'varscan,vardict'.", default="varscan,vardict")
    args = parser.parse_args()

    threads = args.threads or run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()

    if not args.name:
        args.name = guess_sample_name(args.input_fastqs)
        if not args.name:
            sys.exit("Ambiguous sample name")

    if " " in args.name:
        args.name - args.name.replace(" ", "_")

    if args.min_vaf is None:
        args.min_vaf = 0.01 if args.min_family_size == 1 else 0.001

    args.reference = os.path.abspath(args.reference)
    args.input_fastqs = [os.path.abspath(path) for path in args.input_fastqs]
    if args.panel:
        args.panel = os.path.abspath(args.panel)
    if args.vep:
        args.vep = os.path.abspath(args.vep)
    os.chdir(args.output)

    args.reference = (glob.glob(f"{args.reference}/*.fna") + glob.glob(f"{args.reference}/*.fa") + glob.glob(f"{args.reference}/*.fasta") + [args.reference])[0]
    ref_dir = os.path.dirname(args.reference)
    if glob.glob(f"{ref_dir}/*.sa"):
        bwa = "bwa"
    elif glob.glob(f"{ref_dir}/*.0123"):
        bwa = "bwa-mem2"
    else:
        sys.exit("Invalid bwa indexes")
    targets_bedfile = (glob.glob(f"{args.panel}/*.bed") + [None])[0] if args.panel else ""
    stats = f"{args.name}.stats.json"
    pipe = Pipe()


    # Remove umis and do some basic fastq qc
    interleaved_fastq = f"{args.name}.interleaved.fastq"
    command = ["udini", "--output", interleaved_fastq,
                        "--stats", stats,
                        "--umi", args.umi]
    if args.interleaved:
        command.append("--interleaved")
    pipe(command + args.input_fastqs)


    base_sam = f"{args.name}.base.sam"
    with open(base_sam, "wb") as f_out:
        pipe([bwa, "mem", "-t", threads, 
                          "-p", # interleaved paired end fastq
                          "-C", # Append fastq comment to sam
                          args.reference, 
                          interleaved_fastq], stdout=f_out)
    os.unlink(interleaved_fastq)


    sorted_sam = f"{args.name}.sorted.sam"
    pipe(["samtools", "sort", "-o", sorted_sam,
                              "-@", threads,
                              base_sam])
    os.unlink(base_sam)

    if args.sam_only:
        return


    deduplicated_fastq = f"{args.name}.deduplicated.fastq"
    pipe(["elduderino", "--output", deduplicated_fastq,
                        "--stats", stats,
                        "--min-family-size", args.min_family_size,
                        "--umi", args.umi,
                        sorted_sam])
    os.unlink(sorted_sam)


    deduplicated_sam = f"{args.name}.deduplicated.sam"
    with open(deduplicated_sam, "wb") as f_out:
        pipe([bwa, "mem", "-t", threads, 
                          "-p", # interleaved paired end fastq
                          "-C", # Append fastq comment to sam
                          "-Y", # Soft clip non-primary reads
                          args.reference, 
                          deduplicated_fastq], stdout=f_out)
    os.unlink(deduplicated_fastq)


    namesorted_sam = f"{args.name}.namesorted.sam"
    pipe(["samtools", "sort", "-n", # sort by name
                              "-o", namesorted_sam,
                              "-@", threads,
                              deduplicated_sam])
    os.unlink(deduplicated_sam)


    pipe(["size", "--stats", stats,
                  "--rnames", args.sizes,
                  "--output", f"{args.name}.sizes.pdf",
                  namesorted_sam])


    ontarget_sam = f"{args.name}.ontarget.sam"
    pipe(["ontarget", "--output", ontarget_sam,
                      "--bed", targets_bedfile,
                      "--stats", stats,
                      "--cnv", args.cnv,
                      "--threads", threads,
                      namesorted_sam])
    os.unlink(namesorted_sam)


    namesorted_sam = f"{args.name}.namesorted.sam"
    pipe(["samtools", "sort", "-n", # sort by name
                              "-o", namesorted_sam,
                              "-@", threads, 
                              ontarget_sam])
    os.unlink(ontarget_sam)


    fixed_sam = f"{args.name}.fixed.sam"
    pipe(["samtools", "fixmate", namesorted_sam, fixed_sam])
    os.unlink(namesorted_sam)


    if args.translocations:
        pipe(["breakpoint", "--output", f"{args.name}.translocations.tsv",
                            fixed_sam])


    no_read_groups_bam = f"{args.name}.no_read_groups.bam"
    pipe(["samtools", "sort", "-o", no_read_groups_bam,
                              "-@", threads, 
                              fixed_sam])
    os.unlink(fixed_sam)


    bam = f"{args.name}.bam"
    # This step is only required to satisfy Mutect2 and possibly other gatk tools
    pipe(["gatk", "AddOrReplaceReadGroups", f"I={no_read_groups_bam}", f"O={bam}", "LB=lb", "PL=ILLUMINA", "PU=pu", f"SM={args.name}"])
    os.unlink(no_read_groups_bam)


    pipe(["samtools", "index", bam])


    if args.panel:
        pipe(["covermi_stats", "--panel", args.panel,
                               "--output", f"{args.name}.covermi.pdf",
                               "--stats", stats,
                               bam])


    pipe(["call_variants", "--reference", args.reference,
                           "--callers", args.callers,
                           "--name", args.name,
                           "--panel", args.panel,
                           "--vep", args.vep,
                           "--min-vaf", args.min_vaf,
                           "--min-alt-reads", args.min_family_size,
                           "--output", ".", # We have already changed directory into the current directory
                           "--threads", threads])


    #vaf_plot = f"{args.name}.vaf.pdf"
    pipe(["vcf_stats", f"{args.name}.vardict.vcf", # May need to change this depending on variant caller performance
                       "--stats", stats])
                       #"--output", vaf_plot])

    print(pipe.durations, file=sys.stderr, flush=True)



def main():
    try:
        cfpipeline2()
    except OSError as e:
        # File input/output error. This is not an unexpected error therfore
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

