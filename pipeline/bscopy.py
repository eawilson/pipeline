import os
import glob
import argparse
import sys
import re

from pipeline import s3_put, s3_list



def bscopy(bsruns, project, include=None, exclude=None, dry_run=False, rename_run=""):
    if rename_run and len(bsruns) > 1:
        sys.exit("--rename-run can only be specified if a single run is selected")

    for bsrun in bsruns:
        awsrun = rename_run or bsrun
        stem = os.path.expanduser(f"~/basespace/Projects/{bsrun}/Samples")
        if not os.path.exists(stem):
            sys.exit(f"Basespace run {bsrun} does not exist")

        for sample in os.listdir(stem):
            if sample.startswith("."):
                continue
            if (include is not None and not re.search(include, sample)) or (exclude is not None and re.search(exclude, sample)):
                print(f"Skipping {sample}")
                continue
            
            prefix = f"projects/{project}/samples/{awsrun}/{sample}"
            
            uploaded = s3_list("omdc-data", prefix)
            for fastq in glob.glob(f"{stem}/{sample}/Files/*fastq.gz"):
                key = "{}/{}".format(prefix.rstrip("/"), os.path.basename(fastq))
                if key in uploaded:
                    print(f"Already uploaded {fastq}")
                elif not dry_run:
                    s3_put("omdc-data", fastq, prefix=prefix)
                else:
                    print(f"{fastq} -> {prefix}")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bsruns', nargs="+", help="Basespace runs.")
    parser.add_argument("-p", "--project", help="AWS project.", required=True)
    parser.add_argument("-r", "--rename-run", help="Rename run when copying to AWS..", default=argparse.SUPPRESS)
    parser.add_argument('-i', "--include", help="Only include fastqs matching this regular expression.", default=argparse.SUPPRESS)
    parser.add_argument('-e', "--exclude", help="Exclude fastqs matching this regular expression.", default=argparse.SUPPRESS)
    parser.add_argument('-d', "--dry-run", help="Do not actually upload anything to aws.", action="store_true", default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        bscopy(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

