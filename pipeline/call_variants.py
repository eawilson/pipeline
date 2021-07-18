#!/usr/bin/env python3

import os
import pdb
import sys
import argparse
import glob

from pipeline import run, Pipe



def call_variants():
    """Cell free pipeline2 variant calling.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam', help="Path of the input bam file.")
    parser.add_argument("-r", "--reference", help="Path to reference genome or containing directory.", required=True)
    parser.add_argument("-C", "--callers", help="Variant callers to use. Valid values are varscan, vardict and mutect2. Defaults to 'varscan,vardict'.", default="varscan,vardict")
    parser.add_argument("-n", "--name", help="Sample name used to name output files. Will be guessed from input fastq if not provided", default="")
    parser.add_argument("-p", "--panel", help="Path to covermi panel which must contain targets bedfile. Required for annotation.", default="")
    parser.add_argument("-v", "--vep", help="Path to vep cache. Required for annotation.", default="")
    parser.add_argument("-f", "--min-vaf", help="Minimum variant allele frequency for a variant to be called.", type=float, default=0)
    parser.add_argument("-a", "--min-alt-reads", help="Minimum number of alt reads for a variant to be called.", type=int, default=2)
    parser.add_argument("-o", "--output", help="Path to write output files to.", default=".")
    parser.add_argument("-t", "--threads", help="Number of threads to use, defaults to all available threads if not specified.", type=int, default=None)
    args = parser.parse_args()

    threads = args.threads or run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()

    if not args.name:
        args.name = os.path.basename(args.input_bam).split(".")[0]

    args.callers = args.callers.lower().replace(",", " ").split()
    for caller in args.callers:
        if caller not in ("varscan", "vardict", "mutect2"):
            sys.exit(f"{caller} is not a recognised variant caller")

    args.reference = os.path.abspath(args.reference)
    args.input_bam = os.path.abspath(args.input_bam)
    if args.panel:
        args.panel = os.path.abspath(args.panel)
    if args.vep:
        args.vep = os.path.abspath(args.vep)
    os.chdir(args.output)

    args.reference = (glob.glob(f"{args.reference}/*.fna") + glob.glob(f"{args.reference}/*.fa") + glob.glob(f"{args.reference}/*.fasta") + [args.reference])[0]
    pipe = Pipe()

    targets_bedfile = glob.glob(f"{args.panel}/*.bed") if args.panel else []
    targets_bedfile = targets_bedfile[0] if len(targets_bedfile) == 1 else ""

    if "vardict" in args.callers and not targets_bedfile:
        sys.exit(f"No targets bedfile found (required by vardict)")
    if "mutect2" in args.callers and not os.path.exists(f"{args.input_bam}.bai"):
        sys.exit(f"No index found for {args.input_bam} (required by mutect2)")

    ###############################################################################################################
    ### VARSCAN                                                                                                 ###
    ###############################################################################################################
    if "varscan" in args.callers:
        mpileup = f"{args.name}.mpileup"
        pipe(["samtools", "mpileup", "-o", mpileup,
                                     "-f", args.reference,
                                     "-A",
                                     "-B",
                                     "-q", "10",
                                     "-d", "10000000",
                                     args.input_bam])

        pvalue_vcf = f"{args.name}.pvalue.vcf"
        with open(pvalue_vcf, "wb") as f_out:
            pipe(["varscan", "mpileup2cns", mpileup,
                                            "--variants",
                                            "--output-vcf", "1",
                                            "--min-coverage", "1",
                                            "--min-var-freq", args.min_vaf,
                                            "--min-avg-qual", "20",
                                            "--min-reads2", args.min_alt_reads,
                                            "--p-value", "0.05",
                                            "--strand-filter", "1"], stdout=f_out)
        os.unlink(mpileup)

        vcf = f"{args.name}.varscan.unfiltered.vcf" if targets_bedfile else f"{args.name}.varscan.vcf"
        pipe(["postprocess_varscan_vcf", pvalue_vcf, "--output", vcf])
        os.unlink(pvalue_vcf)

        if targets_bedfile:
            unfiltered_vcf = vcf
            vcf = f"{args.name}.varscan.vcf"
            pipe(["filter_vcf", unfiltered_vcf, "--output", vcf, "--bed", targets_bedfile])
            os.unlink(unfiltered_vcf)

        if args.vep and args.panel:
            pipe(["annotate_panel", "--vep", args.vep,
                                    "--output", f"{args.name}.varscan.annotation.tsv",
                                    "--reference", args.reference,
                                    "--threads", threads,
                                    "--panel", args.panel,
                                    vcf])


    ###############################################################################################################
    ### VARDICT                                                                                                 ###
    ###############################################################################################################
    if "vardict" in args.callers:
        vardict_table = f"{args.name}.vardict.tsv"
        with open(vardict_table, "wb") as f_out:
            pipe(["vardictjava", "-K", # include Ns in depth calculation
                                 "-deldupvar", # variants are only called if start position is inside the region interest
                                 "-G", args.reference,
                                 "-N", args.name,
                                 "-b", args.input_bam,
                                 "-Q", "10",
                                 "-f", args.min_vaf,
                                 "-r", args.min_alt_reads,
                                 "-th", threads,
                                 "-u", # count mate pair overlap only once
                                 "-fisher", # perform work of teststrandbias.R
                                 targets_bedfile], stdout=f_out)

        unfiltered_vcf = f"{args.name}.vardict.unfiltered.vcf"
        with open(vardict_table, "rb") as f_in:
            with open(unfiltered_vcf, "wb") as f_out:
                pipe(["var2vcf_valid.pl", "-A", # output all variants at same position
                                          "-f", args.min_vaf,
                                          "-N", args.name], stdin=f_in, stdout=f_out)
        os.unlink(vardict_table)
        
        vcf = f"{args.name}.vardict.vcf"
        # Although vardict take the targets bedfile as an argument is does call occasional variants just outside 
        pipe(["filter_vcf", unfiltered_vcf, "--output", vcf, "--bed", targets_bedfile])
        #os.unlink(unfiltered_vcf)

        if args.vep and args.panel:
            pipe(["annotate_panel", "--vep", args.vep,
                                    "--output", f"{args.name}.vardict.annotation.tsv",
                                    "--reference", args.reference,
                                    "--threads", threads,
                                    "--panel", args.panel,
                                    vcf])


    ###############################################################################################################
    ### MUTECT2                                                                                                 ###
    ###############################################################################################################
    if "mutect2" in args.callers:
        unmutectfiltered_vcf = f"{args.name}.unmutectfiltered.mutect2.vcf"
        pipe(["gatk", "Mutect2", "-R", args.reference,
                                 "-I", args.input_bam,
                                 "-O", unmutectfiltered_vcf,
                                 "--create-output-variant-index", "false",
                                 "--max-reads-per-alignment-start", "0",
                                 "--disable-read-filter", "NotDuplicateReadFilter",
                                 "--disable-read-filter", "GoodCigarReadFilter"])

        multiallelic_vcf = f"{args.name}.multiallelic.mutect2.vcf"
        pipe(["gatk", "FilterMutectCalls", "-R", args.reference,
                                           "-V", unmutectfiltered_vcf,
                                           "-O", multiallelic_vcf,
                                           "--filtering-stats", "false",
                                           "--create-output-variant-index", "false"])
        os.unlink(unmutectfiltered_vcf)
        os.unlink(f"{unmutectfiltered_vcf}.stats")

        vcf = f"{args.name}.mutect2.unfiltered.vcf" if targets_bedfile else f"{args.name}.mutect2.vcf"
        pipe(["postprocess_mutect2_vcf", "--output", vcf,
                                         "--min-alt-reads", args.min_alt_reads,
                                         "--min-vaf", args.min_vaf,
                                         multiallelic_vcf])
        os.unlink(multiallelic_vcf)

        if targets_bedfile:
            unfiltered_vcf = vcf
            vcf = f"{args.name}.mutect2.vcf"
            pipe(["filter_vcf", unfiltered_vcf, "--output", vcf, "--bed", targets_bedfile])
            os.unlink(unfiltered_vcf)

        if args.vep and args.panel:
            pipe(["annotate_panel", "--vep", args.vep,
                                    "--output", f"{args.name}.mutect2.annotation.tsv",
                                    "--reference", args.reference,
                                    "--threads", threads,
                                    "--panel", args.panel,
                                    vcf])


    print(pipe.durations, file=sys.stderr, flush=True)



def main():
    try:
        call_variants()
    except OSError as e:
        # File input/output error. This is not an unexpected error therfore
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

