import os
import glob
import pdb
from pipeline import run, illumina_readgroup, pipe, create_report, command_line_arguments, sample_name, vcf_pvalue_2_phred
from covermi import Panel, covermimain



def cfpipeline(fastq_r1, fastq_r2, panel, genome_dir=".", dest_dir="."):
    sample = sample_name(fastq_r1, fastq_r2)
    threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    prop = Panel(panel).properties
    reference_fasta = glob.glob(os.path.join(genome_dir, prop["assembly"], "*", "*.fna"))[0]
    
    dude_options = ["--allowed", "3",
                    "--debug_dump_umis",
                    "--min_family_size", prop.get("min_family_size", "1")]
    if prop.get("thruplex", False):
        dude_options += ["--thruplex"]
    pipe(["dude", fastq_r1, fastq_r2] + dude_options)
    deduped_fastqs = ["{}.deduped.fastq".format(fastq[:-6]) for fastq in (fastq_r1, fastq_r2)]
    pipe(["fastuniq", "-i", fastq_r1, fastq_r2, "-o", deduped_fastqs[0], "-p", deduped_fastqs[1]])
    
    sam_file = "{}.sam".format(sample)
    with open(sam_file, "wb") as f_out:
        pipe(["bwa", "mem", "-t", threads, "-R", illumina_readgroup(deduped_fastqs[0]), reference_fasta] + deduped_fastqs, stdout=f_out)
    for fastq in deduped_fastqs:
        os.unlink(fastq)

    unsorted_bam_file = "{}.unsorted.bam".format(sample)
    pipe(["samtools", "fixmate", "-O", "bam", sam_file, unsorted_bam_file])
    os.unlink(sam_file)

    bam_file = "{}.bam".format(sample)
    pipe(["samtools", "sort", "-O", "bam", "-o", bam_file, "-T", "temp", "-@", threads, unsorted_bam_file])
    os.unlink(unsorted_bam_file)

    pipe(["samtools", "index", bam_file])

    mpileup_file = "{}.pileup".format(sample)
    mpileup_options = ["-A",
                       #"-B",
                       "-q", "10",
                       "-d", "10000000"]
    pipe(["samtools", "mpileup", "-o", mpileup_file, "-f", reference_fasta] + mpileup_options + [bam_file])
    
    pvalue_vcf_file = "{}.pvalue.vcf".format(sample)
    varscan_options = ["--min-coverage", "1",
                       "--min-var-freq", "0", 
                       "--min-avg-qual", "20",
                       "--min-reads2", "1",
                       "--p-value", "0.001"]
                    #"--min-coverage-normal", "1",  
                    #"--min-coverage-tumor", "1",  
                    #"--min-freq-for-hom","0.75", 
                    #"--somatic-p-value", "0.05",
                    #"--strand-filter", "1",
                    #"--validation", "1"

    with open(temp_vcf_file, "wb") as f_out:
        pipe(["varscan", "mpileup2cns", mpileup_file, "--variants", "--output-vcf", "1"] + varscan_options, stdout=f_out) 
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
    pipe(["perl", vep_script, "-i", vcf_file, "-o", vepjson_file, "--fasta", reference_fasta] + vep_options)
    annotation_file = create_report(vepjson_file, panel)
    os.unlink(vepjson_file)

    covermi_dir = covermimain(panel, "", bam_path=bam_file)
    covermi_file = "{}.covermi.tar.gz".format(sample)
    run(["tar", "cvzf", covermi_file, covermi_dir])
    run(["rm", "-r", covermi_dir])



def main():
    cfpipeline(**command_line_arguments())



if __name__ == "__main__":
    main()

