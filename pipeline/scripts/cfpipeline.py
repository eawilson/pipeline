import os
import pdb
import sys
import argparse
from pipeline import run, illumina_readgroup, pipe, create_report, \
        command_line_arguments, sample_name, vcf_pvalue_2_phred, \
        ungzip_and_combine_illumina_fastqs, fasta_path, \
        sam_remove_offtarget
from covermi import Panel, covermimain



def cfpipeline(fastqs, panel, genome, min_family_size):
    panel_name = panel

    threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    sample = sample_name(fastqs)
    merged_fastqs = ungzip_and_combine_illumina_fastqs(fastqs, sample)
    if merged_fastqs is None:
        merged_fastqs = fastqs
        fastqs_to_be_deleted = ()
    else:
        fastqs_to_be_deleted = merged_fastqs
    
    panel = Panel(panel_name)
    prop = panel.properties
    reference_fasta = fasta_path(genome)
    
    dude_options = ["--allowed", "3",
                    "--min_family_size", min_family_size]
    if prop.get("thruplex", False):
        dude_options += ["--thruplex"]
    pipe(["dude"] + merged_fastqs + dude_options)
    deduped_fastqs = ["{}.deduped.fastq".format(fastq[:-6])
                        for fastq in merged_fastqs]
    for fastq in fastqs_to_be_deleted:
        os.unlink(fastq)
    #pipe(["fastuniq", "-i", fastq_r1, fastq_r2, "-o", deduped_fastqs[0], "-p", deduped_fastqs[1]])
    
    sam_file = "{}.sam".format(sample)
    with open(sam_file, "wb") as f_out:
        pipe(["bwa", "mem", "-t", threads, 
                            "-R", illumina_readgroup(deduped_fastqs[0]),
                            reference_fasta] + deduped_fastqs, stdout=f_out)
    for fastq in deduped_fastqs:
        os.unlink(fastq)
        
    if "amplicons" in panel:
        ontarget_sam_file = "{}.ontarget.sam".format(sample)
        sam_remove_offtarget(sam_file, ontarget_sam_file, panel.amplicons)
        os.unlink(sam_file)
        sam_file = ontarget_sam_file

    unsorted_bam_file = "{}.unsorted.bam".format(sample)
    pipe(["samtools", "fixmate", "-O", "bam", sam_file, unsorted_bam_file])
    #os.unlink(sam_file)

    bam_file = "{}.bam".format(sample)
    pipe(["samtools", "sort", "-O", "bam",
                              "-o", bam_file,
                              "-T", "temp",
                              "-@", threads, 
                              unsorted_bam_file])
    os.unlink(unsorted_bam_file)

    pipe(["samtools", "index", bam_file])

    mpileup_file = "{}.pileup".format(sample)
    mpileup_options = ["-A",
                       "-B",
                       "-q", "10",
                       "-d", "10000000"]
    pipe(["samtools", "mpileup", "-o", mpileup_file,
                                 "-f", reference_fasta] + 
                                 mpileup_options + [bam_file])
    
    pvalue_vcf_file = "{}.pvalue.vcf".format(sample)
    varscan_options = ["--min-coverage", "1",
                       "--min-var-freq", "0", 
                       "--min-avg-qual", "20",
                       "--min-reads2", "3",
                       "--p-value", "0.05"]
                    #"--min-coverage-normal", "1",
                    #"--min-coverage-tumor", "1",
                    #"--min-freq-for-hom","0.75",
                    #"--somatic-p-value", "0.05",
                    #"--strand-filter", "1",
                    #"--validation", "1"

    with open(pvalue_vcf_file, "wb") as f_out:
        pipe(["varscan", "mpileup2cns", mpileup_file, "--variants", 
                                                      "--output-vcf", "1"] + 
                                                      varscan_options, stdout=f_out)
    os.unlink(mpileup_file)
    vcf_file = "{}.vcf".format(sample)
    vcf_pvalue_2_phred(pvalue_vcf_file, vcf_file)
    os.unlink(pvalue_vcf_file)
    
    vepjson_file = "{}.vep".format(sample)
    vep_options = ["--verbose",
                   "--no_stats",
                   "--format", "vcf",
                   "--fork", threads,
                   "--json",
                   "--offline",
                   "--everything",
                   "--force_overwrite"]
    if prop.get("transcript_source", "refseq") == "refseq":
        vep_options += ["--refseq"]
    if "assembly" in prop:
        vep_options += ["--assembly", prop["assembly"]]
    vep_script = os.path.join(os.path.expanduser("~"), "ensembl-vep/vep")
    pipe(["perl", vep_script, "-i", vcf_file,
                              "-o", vepjson_file,
                              "--fasta", reference_fasta] + vep_options)
    annotation_file = create_report(vepjson_file, panel)
    os.unlink(vepjson_file)

    covermi_dir = covermimain(panel_name, "", bam_path=bam_file)
    covermi_file = "{}.covermi.tar.gz".format(sample)
    run(["tar", "cvzf", covermi_file, covermi_dir])
    run(["rm", "-r", covermi_dir])



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastqs', nargs="+", help="Fastq files.")
    parser.add_argument("-p", "--panel", help="Directory containing panel data.")
    parser.add_argument("-g", "--genome", help="Reference genome.")
    #parser.add_argument("-b", "--bam", help="Matched normal bam to perform tumour/normal subtraction.", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--min-family-size", help="Minimum family size.", dest="min_family_size", default="1")
    args = parser.parse_args()
    cfpipeline(**vars(args))



if __name__ == "__main__":
    main()

