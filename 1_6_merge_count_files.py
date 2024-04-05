"""
This script parsers the 1_5_bam_to_counts.py output files and creates the input formats for DESeq2.

    Args:
        -S, --suffix (str): A recognizable suffix included in the file names to identify the count files to be merged.
        -G, --groups (str): A comma-separated list of the codes in the file names identifying the experiment groups. Groups' codes must be input as a comma-separated string.
        -T, --read_type (str): The type of sequencing, either single-read or pair-end.
        -F, --folder_path (str): The path to the folder where the output files are stored. If folder_path is not provided, it defaults to "06_output/"
        -D, --colData (str): When specified, it creates the colData file for DESeq2.
        -R, --run (str): The control variable for the function to be run (1 to run).
"""

import argparse
from functions.libs import merge_count_files, create_colData, mkdir

# Gets the command line arguments with argparse.
parser = argparse.ArgumentParser()
parser.add_argument("-S", "--suffix")
parser.add_argument("-G", "--groups")
parser.add_argument("-T", "--read_type")
parser.add_argument("-F", "--folder_path", type=str, default="06_output/")
parser.add_argument("-D", "--colData", action="store_true")
parser.add_argument("-R", "--run", type=bool, default=False)
args = vars(parser.parse_args())

# Assign the command line arguments to variables.
suffix, groups, read_type, folder_path, colData, run = (
    args["suffix"],
    args["groups"],
    args["read_type"],
    args["folder_path"],
    args["colData"],
    args["run"],
)

# groups should be input as a comma-separated string.
groups = tuple(map(str, groups.split(",")))

# The number indicates the script that has generated the file (except for 00_log).
mkdir("06_output/")

# Merges all the sample counts to create the count matrix.
merge_count_files(suffix=suffix, run=run, folder_path=folder_path)

# Creates the colData file on request.
if colData:

    create_colData(
        groups=groups,
        suffix=suffix,
        read_type=read_type,
        run=run,
        folder_path=folder_path,
    )
