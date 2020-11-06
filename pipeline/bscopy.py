import os
import glob
import argparse

from pipeline import s3_put, s3_list



def bscopy(bsruns, project, dry_run=False):
    for bsrun in bsruns:
        stem = os.path.expanduser(f"~/basespace/Projects/{bsrun}/Samples")
        for sample in os.listdir(stem):
            if sample.startswith("."):
                continue
            
            #if not sample.endswith("Combined"):
                #continue
            
    #        real_sample = sample[:-8].upper()

            prefix = f"projects/{project}/samples/{bsrun}"
            if True:#not s3_list("omdc-data", prefix):
                for fastq in glob.glob(f"{stem}/{sample}/Files/*fastq.gz"):
                    if not dry_run:
                        s3_put("omdc-data", fastq, prefix=prefix)
                    elsw:
                        print(fastq)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bsruns', nargs="+", help="Basespace runs.")
    parser.add_argument("-p", "--project", help="AWS project.", required=True)
    parser.add_argument('-d', "--dry-run", action="store_const", const=True, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        bscopy(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

