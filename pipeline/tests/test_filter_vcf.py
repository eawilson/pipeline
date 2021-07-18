import subprocess
import os



source_dir = os.path.dirname(__file__)
bed = os.path.join(source_dir, "filter_vcf.bed")
vcf = os.path.join(source_dir, "filter_vcf.vcf")

with open(vcf, "rt") as f_in:
    overlaps = f_in.read().count("overlap")

cmd = ["filter_vcf", vcf, "--output", "-", "--bed", bed]
cp = subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True)

if cp.stdout.count("overlap") == overlaps:
    print("passed")
else:
    print("failed")





