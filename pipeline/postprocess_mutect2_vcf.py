import math
import argparse
import sys
from math import log
import pdb


ALT = 4
QUAL = 5
FILTER = 6
INFO = 7
FORMAT_KEYS = 8
FORMAT_VALS = 9



def tlod2phred(tlod):
    # TLOD="Log 10 likelihood ratio score of variant existing versus not existing"
    try:
        likelyhood_ration= 10 ** float(tlod)
    except OverflowError:
        return 93
    phred = -10 * log(1 / (1 + likelyhood_ration), 10)
    return str(min(int(phred), 93))



def postprocess_mutect2_vcf(input_vcf, output_vcf=None, min_vaf=0, min_alt_reads=0):
    """ Mutect2 produces vcfs with a number of feature that are not
        compatible with downstream process. This script will output
        a new vcf with these corrected.
        1) Split multi-allelic lines into one line per alt call.
        2) Remove multiallelic from the filters column and replace
            with PASS.
        3) Calculate a phred quality score from the tlod score. ?Is
            this correct as the qualities are all very high.
        Filter output by min-alt-reads and min-vaf if provided.
        NOTE - Will not work for vcf with multiple samples.
    """

    if not input_vcf.endswith(".vcf"):
        sys.exit("Input must be a .vcf file")
    if output_vcf is None:
        output_vcf = "{}.postprocessed.vcf".format(input_vcf[:-4])
    elif not output_vcf.endswith(".vcf"):
        sys.exit("Output must be a .vcf file")

    sources = []
    info_numbers = {}
    format_numbers = {}

    with open(input_vcf, "rt") as f_in:
        with open(output_vcf, "wt") as f_out:
            for row in f_in:
                row = row.split("\t")
                if row[0].startswith("#"):
                    if row[0].startswith("##source="):
                        sources.append(row[0][9:])

                    if row[0].startswith("##INFO="):
                        info_numbers[row[0].split("ID=")[1].split(",")[0]] = row[0].split("Number=")[1].split(",")[0]

                    if row[0].startswith("##FORMAT="):
                        format_numbers[row[0].split("ID=")[1].split(",")[0]] = row[0].split("Number=")[1].split(",")[0]

                    f_out.write("\t".join(row))
                    continue
                
                if len(row) > 10:
                    sys.exit("Multi-sample vcf files not supported")

                alts = row[ALT].split(",")


                quals = [tlod2phred(tlod) for tlod in row[INFO].split("TLOD=")[1].split(";")[0].split(",")]
                if len(quals) != len(alts):
                    sys.exit("Invalid vcf, mismatch between number of alleles and number of quality scores")


                default = row[FILTER]
                if len(alts) > 1:
                    default = ";".join(val for val in default.split(";") if val != "multiallelic")
                    if default == "":
                        default = "PASS"
                filters = [(fil if fil != "SITE" else default) for fil in row[INFO].split("AS_FilterStatus=")[1].split(";")[0].split("|")]



                if len(alts) == 1:
                    infos = [row[INFO]]

                else:
                    infos = [[] for alt in alts]
                    for keyval in row[INFO].split(";"):
                        key = keyval.split("=")[0] if "=" in keyval else keyval

                        if info_numbers[key] not in ("R", "A"):
                            for i, alt in enumerate(alts):
                                infos[i].append(keyval)
                            continue

                        key, val = keyval.split("=")
                        sep = "|"
                        vals = val.split(sep)
                        if len(vals) != len(alts) + bool(info_numbers[key] == "R"):
                            sep = ","
                            vals = val.split(sep)
                            if len(vals) != len(alts) + bool(info_numbers[key] == "R"):
                                sys.exit(f"Invalid vcf, incorrect number of values for INFO {key}")

                        if info_numbers[key] == "R":
                            key = "{}={}{}".format(key, vals[0], sep)
                            vals = vals[1:]
                        else:
                            key = f"{key}="

                        for i, val in enumerate(vals):
                            infos[i].append(key + val)

                    infos = [";".join(keyvals) for keyvals in infos]


                allelic_depths = None
                allelic_fractions = None
                formats = [[] for alt in alts]
                for key, val in zip(row[FORMAT_KEYS].split(":"), row[FORMAT_VALS].split(":")):

                    if format_numbers[key] not in ("R", "A"):
                        for i, alt in enumerate(alts):
                            formats[i].append(val)
                        continue

                    vals = val.split(",")
                    if len(vals) != len(alts) + bool(format_numbers[key] == "R"):
                        sys.exit(f"Invalid vcf, incorrect number of values for FORMAT {key}")

                    if format_numbers[key] == "R":
                        prefix = "{},".format(vals[0])
                        vals = vals[1:]
                    else:
                        prefix = ""

                    if key == "AD":
                        allelic_depths = [int(val) for val in vals]

                    if key == "AF":
                        allelic_fractions = [float(val) for val in vals]

                    for i, val in enumerate(vals):
                        formats[i].append(prefix + val)

                formats = [":".join(vals) for vals in formats]



                for alt, qual, filt, info, fmat, alt_reads, vaf in zip(alts, quals, filters, infos, formats, allelic_depths, allelic_fractions):
                    if min_alt_reads and alt_reads < min_alt_reads:
                        continue
                    if min_vaf and vaf < min_vaf:
                        continue

                    row[ALT] = alt
                    row[QUAL] = qual
                    row[FILTER] = filt
                    row[INFO] = info
                    row[FORMAT_VALS] = fmat
                    f_out.write("\t".join(row))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_vcf', help="Input vcf file.")
    parser.add_argument("-o", "--output", help="Output vcf file.", dest="output_vcf", default=argparse.SUPPRESS)
    parser.add_argument("-m", "--min-vaf", help="All calls with a vaf less than this will be filtered.", type=float, default=argparse.SUPPRESS)
    parser.add_argument("-a", "--min-alt-reads", help="All calls with a number of alt reads less than this will be filtered.", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        postprocess_mutect2_vcf(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

