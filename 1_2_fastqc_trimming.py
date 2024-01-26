import argparse
import json
from functions.libs import (
    list_dir_files,
    get_sample_name,
    mkdir,
    eval_fastq_files,
    trimming_files,
    get_stats_fastq_files,
    trimming_files_slow,
)

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--input-dir")
parser.add_argument("-A", "--adapter")
parser.add_argument("-R", "--run")
parser.add_argument(
    "-L", "--slow"
)  # Much slower processing, but less memory-intensive.
parser.add_argument("-P", "--processes")
args = vars(parser.parse_args())

input_dir, adapter, run, slow, processes = (
    args["input_dir"],
    args["adapter"],
    args["run"],
    args["slow"],
    args["processes"],
)

if processes:
    processes = int(processes)

filenames = list_dir_files(input_dir, "fastq.gz")
sample_names = get_sample_name(filenames)

sample_dict = {}
for sample_name in sample_names:
    fastq_file_r1 = [x for x in filenames if sample_name in x and "_R1_" in x][0]
    sample_dict[sample_name] = fastq_file_r1

mkdir("FastQC")
mkdir("FastQC/Raw")
mkdir("FastQC/Trim")
mkdir("00_log")
mkdir("02_trim")
mkdir("05_plot")


if slow == "1":
    sample_dict = trimming_files_slow(sample_dict, adapter, run)
else:
    eval_fastq_files(sample_dict, "FastQC/Raw", adapter, run, processes)
    sample_dict = trimming_files(sample_dict, adapter, run)
    eval_fastq_files(sample_dict, "FastQC/Trim", "None", run, processes)
    get_stats_fastq_files(sample_dict, run, processes)

with open("00_log/1_2_fastq.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)
