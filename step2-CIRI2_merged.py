# Version 1.1
# Date: 10/04/2017

import argparse
import sys
import os
import glob
import pandas as pd
import subprocess

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates circRNA predict scripts:CIRI .")
parser.add_argument("-i", "--inputDiretory", help="Directory with subdirectories include *.CIRI.gtf (force).")
parser.add_argument("-o", "--outputDirectory", help="Output directory with merged circRNA results (force).")
parser.add_argument("-p", "--prefix", help="output file's prefix.[default:combinded_circRNA].",
                    default="combinded_circRNA")
args = parser.parse_args()

# Process the command line arguments.

inputDiretory = os.path.abspath(args.inputDiretory)
prefix = args.prefix
outputDirectory = os.path.abspath(args.outputDirectory)

# Creat outputDirectory if it is not exist
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
os.chdir(outputDirectory)


# define varible
# global merge

# define function for
def circRNA_merged(input, output, name):
    files = glob.glob(input + "/*/*.CIRI.gtf")
    files.sort()
    if os.path.exists(os.path.join(output, name + ".bed")):
        os.remove(os.path.join(output, name + ".bed"))
    for i in range(0, len(files), 1):
        filename = os.path.basename(files[i]).replace(".CIRI.gtf", "")
        odata = pd.read_csv(files[i], header=0, sep='\t')
        for j in range(0, len(odata)):
            odata.iloc[j, 2] = odata.iloc[j, 2] - 1
        bed = pd.DataFrame({'0_chr': odata.iloc[:, 1],
                            '1_start': odata.iloc[:, 2],
                            '2_end': odata.iloc[:, 3],
                            '3_circRNA_id': odata.iloc[:, 0],
                            '4_junction_reads': odata.iloc[:, 4],
                            '5_strand': odata.iloc[:, 10]})

        if not os.path.exists(os.path.join(output, filename, filename + ".CIRI.bed")):
            bed.to_csv(os.path.join(output, filename, filename + ".CIRI.bed"), sep='\t', header=False, encoding='utf-8',
                       index=False)
        for k in range(0, len(bed)):
            bed.iloc[k, 4] = "0"
        if i == 0:
            merge = bed
        if i > 0:
            merge = pd.concat([merge, bed], axis=0, ignore_index=True)

    merge.drop_duplicates(keep='first', subset='3_circRNA_id', inplace=True)
    merge.sort_values(['0_chr', '1_start', '2_end'], ascending=[True, True, True], inplace=True)
    merge.reset_index(drop=True, inplace=True)
    number = len(str(len(merge)))

    for l in range(0, len(merge)):
        new_name = str(str(l + 1).zfill(number))
        merge.iloc[l, 3] = str(
            str("circRNA_") + new_name + "|" + str(merge.iloc[l, 0]) + ":" + str(merge.iloc[l, 1] + 1) + "_" + str(
                merge.iloc[l, 2]) + "_" + str(merge.iloc[l, 5]))
    merge.to_csv(os.path.join(output, name + ".bed"), sep='\t', header=False, encoding='utf-8', index=False)

    files = glob.glob(input + "/*/*.CIRI.bed")
    files.sort()
    if os.path.exists(os.path.join(output, "counts_junction_reads.xls")):
        os.remove(os.path.join(output, "counts_junction_reads.xls"))

    combind_bed = pd.read_csv(os.path.join(output, name + ".bed"), sep='\t', names=["chr", "start", "end", "id", "number", "strand"])
    combind_bed['IDS'] = 'IDS'
    for l in range(0, len(combind_bed)):
        combind_bed.iloc[l, 6] = str(
            str(combind_bed.iloc[l, 0]) + ":" + str(combind_bed.iloc[l, 1] + 1) + "|" + str(combind_bed.iloc[l, 2]))

    for i in range(0, len(files), 1):
        filename = os.path.basename(files[i]).replace(".CIRI.bed", "").replace("Sample_", "")
        odata = pd.read_csv(files[i], sep='\t', index_col=3, names=["chr", "start", "end", filename, "strand"])
        odata.drop(["chr", "start", "end", "strand"], axis=1, inplace=True)
        if i == 0:
            merge_reads = combind_bed.join(odata, on=['IDS']).fillna("0")
        if i >= 1:
            merge_reads = merge_reads.join(odata, on=['IDS']).fillna("0")
    merge_reads.reset_index(drop=True, inplace=True)
    merge_reads.drop(["chr", "start", "end", "IDS", "number", "strand"], axis=1, inplace=True)
    merge_reads.replace("Sample_", "")
    merge_reads.to_csv(os.path.join(output, "counts_junction_reads.xls"), sep='\t', header=True, encoding='utf-8',
                       index=False)


if __name__ == "__main__":
    # Bwa +CIRI execution
    circRNA_merged(inputDiretory, outputDirectory, prefix)
