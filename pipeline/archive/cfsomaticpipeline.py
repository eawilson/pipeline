import os
import pdb
import sys
import argparse
from pipeline import run, illumina_readgroup, pipe, create_report, \
        command_line_arguments, sample_name, vcf_pvalue_2_phred, \
        ungzip_and_combine_illumina_fastqs, fasta_path, \
        sam_remove_offtarget
from covermi import Panel



def cfsomaticpipeline(bam_files, panel, genome):
    if len(bam_files) != 2:
        raise RuntimeError("Must be exactly two bams.")
    
    panel_name = panel

    threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    sample = bam_files[1][:-4] +".somatic"
    panel = Panel(panel_name)
    prop = panel.properties
    reference_fasta = fasta_path(genome)

    mpileup_file = "{}.pileup".format(sample)
    mpileup_options = ["-A",
                    "-B",
                    "-q", "10",
                    "-d", "10000000"]
    pipe(["samtools", "mpileup", "-o", mpileup_file,
                                "-f", reference_fasta] + 
                                mpileup_options + bam_files)
    
    pvalue_vcf_file = "{}.pvalue.vcf".format(sample)
    indel_file = "{}.pvalue.indel.vcf".format(sample)
    varscan_options = ["--output-snp", pvalue_vcf_file,
                       "--output-indel", indel_file,
                       "--mpileup", "1",
                       "--output-vcf", "1",
                       "--min-coverage", "1",
                       "--min-var-freq", "0", 
                       "--min-avg-qual", "20",
                       "--min-reads2", "3",
                       "--p-value", "0.05",
                       "--min-coverage-normal", "1",
                       "--min-coverage-tumor", "1",
                       "--min-freq-for-hom","0.75",
                       "--somatic-p-value", "0.05",
                       "--strand-filter", "1"]
    pipe(["varscan", "somatic", mpileup_file] + varscan_options)
    os.unlink(mpileup_file)

    with open(pvalue_vcf_file, "a+t") as f_out:
        with open(indel_file, "rt") as f_in:
            for row in f_in:
                if not row.startswith("#"):
                    f_out.write(row)
    os.unlink(indel_file)
        
    vcf_file = "{}.vcf".format(sample)
    vcf_pvalue_2_phred(pvalue_vcf_file, vcf_file)
    os.unlink(pvalue_vcf_file)
    
    vepjson_file = "{}.vep".format(sample)
    vep_options = ["--no_stats",
                   "--format", "vcf",
                   "--fork", threads,
                   "--json",
                   "--offline",
                   "--everything",
                   "--warning_file", "STDERR",
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



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_files', nargs="+", help="Normal bam followed by tumour bam.")
    parser.add_argument("-p", "--panel", help="Directory containing panel data.")
    parser.add_argument("-g", "--genome", help="Reference genome.")
    args = parser.parse_args()
    cfsomaticpipeline(**vars(args))



if __name__ == "__main__":
    main()

