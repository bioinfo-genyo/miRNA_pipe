import argparse
from functions.libs import merge_count_files, create_colData

# This function parsers the 1_5_bam_to_counts.py output files and creates the input formats for DESeq2.
# Suffix
# If colData is set to 1, it also creates the colData file.
# Codes should be input as a comma-separated string, being recognizable patterns in the file names to identify the sample group.
# If folder_path is not provided, it defaults to "04_counts/"
# read_type is the type of sequencing, either single-read or pair-end.
# run is a control for the function to be run.

parser = argparse.ArgumentParser()
parser.add_argument("-S", "--suffix")
parser.add_argument("-C", "--codes")
parser.add_argument("-T", "--read_type")
parser.add_argument("-F", "--folder_path")
parser.add_argument("-D", "--colData")
parser.add_argument("-R", "--run")
args = vars(parser.parse_args())

suffix, codes, read_type, folder_path, colData, run = (
    args["suffix"],
    args["codes"],
    args["read_type"],
    args["folder_path"],
    args["colData"],
    args["run"],
)

# Codes should be input as a comma-separated string.
codes = tuple(map(str, codes.split(",")))

if not folder_path:
    folder_path = "04_counts/"

# Merges all the sample counts to create the count matrix.
merge_count_files(suffix=suffix, run=run, folder_path=folder_path)

# Creates the colData file on request.
if colData == "1":

    create_colData(
        codes=codes,
        suffix=suffix,
        read_type=read_type,
        run=run,
        folder_path=folder_path,
    )
