import os
import pdb
import sys
import argparse
import glob
from pipeline import run, pipe, vcf_pvalue_2_phred, create_report
from covermi import Panel, covermimain



def cfpipeline(sample, input_fastqs, reference, panel, umi=None, vep=None):
    """Cell free pipeline.

    Args:
        sample (str): Sample name, used to name output files.
        input_fastqs (list of str): Paths of input paired fastqs or fastq.gzs,
        reference (str): Path to reference fasta.
        targets_bedfile (str): Path to bedfile of baits/amplicons/regions of interest.
        
    Returns:
        None.
        
    Raises:
    """
    threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    
    reference = (glob.glob(f"{reference}/*.fna") + [reference])[0]
    targets_bedfile = (glob.glob(f"{panel}/*.bed") + [panel])[0]
    
    deduped_interleaved_fastq = f"{sample}.interleaved.fastq"
    stats = f"{sample}.stats.json"
    dude_options = ["--output", deduped_interleaved_fastq,
                    "--stats", stats]
    if umi is not None:
        dude_options += ["--umi", umi]
    pipe(["dude"] + input_fastqs + dude_options)


    unsorted_unfixed_sam = f"{sample}.unsorted.unfixed.sam"
    with open(unsorted_unfixed_sam, "wb") as f_out:
        pipe(["bwa", "mem", "-t", threads, 
                            "-p", 
                            "-C", 
                            reference, 
                            deduped_interleaved_fastq], stdout=f_out)
    os.unlink(deduped_interleaved_fastq)

    namesorted_unfixed_sam = f"{sample}.namesorted.unfixed.sam"
    pipe(["samtools", "sort", "-n",
                              "-o", namesorted_unfixed_sam,
                              "-T", "temp",
                              "-@", threads, 
                              unsorted_unfixed_sam])
    os.unlink(unsorted_unfixed_sam)
    
    namesorted_sam = f"{sample}.namesorted.sam"
    pipe(["samtools", "fixmate", namesorted_unfixed_sam, namesorted_sam])
    os.unlink(namesorted_unfixed_sam)
    
    unfiltered_sam = f"{sample}.unfiltered.sam"
    pipe(["samtools", "sort", "-o", unfiltered_sam,
                              "-T", "temp",
                              "-@", threads, 
                              namesorted_sam])
    os.unlink(namesorted_sam)
    
    sam = f"{sample}.sam"    
    filter_options = ["--output", sam,
                      "--stats", stats,
                      "--min-family-size", "2"]
    if targets_bedfile is not None:
        filter_options += ["--targets", targets_bedfile]
    pipe(["filter_sam", unfiltered_sam] + filter_options)
    os.unlink(unfiltered_sam)

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
    vep_script = os.path.join(os.path.expanduser("~"), "ensembl-vep/vep")
    pipe(["perl", vep_script, "-i", vcf,
                              "-o", vepjson,
                              "--fasta", reference] + vep_options)
    create_report(vepjson, panel)
    os.unlink(vepjson)

    covermi_dir = covermimain(panel, "", bam_path=bam)
    covermi_file = "{}.covermi.tar.gz".format(sample)
    run(["tar", "cvzf", covermi_file, covermi_dir])
    run(["rm", "-r", covermi_dir])



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastqs', nargs="+", help="Fastq files.")
    parser.add_argument("-s", "--sample", help="Sample name.", required=True)
    parser.add_argument("-r", "--reference", help="Reference genome.", required=True)
    parser.add_argument("-p", "--panel", help="Directory containing panel data.", required=True)
    parser.add_argument("-v", "--vep", help="Directory containing vep data.", required=True)
    parser.add_argument("-u", "--umi", help="Umi.", default=argparse.SUPPRESS)
    args = parser.parse_args()
    cfpipeline(**vars(args))



if __name__ == "__main__":
    main()

