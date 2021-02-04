import pdb
import argparse
import csv
import sys


def rename_bed(input, output="output.bed"):
    with open(input, "rt") as f_in:
        reader = csv.reader(f_in, delimiter="\t")
        with open(output, "wt") as f_out:
            writer = csv.writer(f_out, delimiter="\t")
            
            for row in reader:
                row = row[:4]
                row[3] = row[3].split("(")[-1].split(")")[0]
                writer.writerow(row)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help="Input bed file.")
    parser.add_argument("-o", "--output", help="Output bed file.", default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        rename_bed(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

