import os
import subprocess
import pdb


os.chdir("/media/psf/Home/Downloads/han")
fastq_dir = "/media/psf/Home/Downloads/han_fastqs"


all_fastqs = os.listdir(fastq_dir)
samples = list(set(fastq.split("_")[0] for fastq in all_fastqs))

for sample in samples:
    print(sample)
    fastqs = [os.path.join(fastq_dir, fastq) for fastq in all_fastqs if fastq.startswith(sample)]

   # os.mkdir("temp")
    os.chdir("temp")

    with open("{}.cfpipeline.txt".format(sample), "wt") as f:
        subprocess.run(["cfpipeline"] + fastqs +
                        ["--panel", "/home/ed/Data/panels/Head_and_Neck",
                         "--genome", "/home/ed/Data/genomes/GRCh37/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.fna",
                         "--min-family-size", "3"],
                        stderr=f, check=True)
                        
    os.chdir("..")
    os.rename("temp", sample)
    break








