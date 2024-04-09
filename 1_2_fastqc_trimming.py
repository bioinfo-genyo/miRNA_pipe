""" 
This is the first step of the data pre-processing, delete the adapter and UMIs from the reads and trim by quality scores.

    Args:
        -I, --input_dir (str): Input directory.
        -A, --adapter (str): Adapter sequence to remove. Default is the Illumina universal adapter.
        -R, --run (str): Run control variable (1 to run).
        -L, --slow (str): Slow mode. It only uses one thread. Use in case the memory use overwhelms the system capabilities.
        -P, --processes (str): The number of cpu threads to use. Default is 4. If 0 is specified, use the number of samples to maximize parallelization.
        -a, --append_sample_dict (str): If specified, appends an existing sample dictionary to the new one.
"""

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

# Gets the command line arguments with argparse.
parser = argparse.ArgumentParser()
parser.add_argument("-I", "--input_dir", type=str)
parser.add_argument("-A", "--adapter", type=str, default="AGATCGGAAGAG")
parser.add_argument("-R", "--run", type=bool, default=False)
parser.add_argument(
    "-L", "--slow", type=bool, default=False
)  # Much slower processing, but less memory-intensive.
parser.add_argument("-P", "--processes", type=int, default=4)
parser.add_argument("-a", "--append_sample_dict", action="store_true")
args = vars(parser.parse_args())

# Assign the command line arguments to variables.
input_dir, adapter, run, slow, processes, append = (
    args["input_dir"],
    args["adapter"],
    args["run"],
    args["slow"],
    args["processes"],
    args["append_sample_dict"],
)

# Build sample dict. Key is sample name, value is fastq file path. This is what we use to localize the appropiate files for each step.
# When a processing step is performed over the samples, their file name is changed to indicate the process performed.
# The sample dict is stored as a json file that is updated on each step to keep track of the transformed samples' fastq files.
filenames = list_dir_files(input_dir, "fastq.gz")
sample_names = get_sample_name(filenames)

sample_dict = {}
for sample_name in sample_names:
    fastq_file_r1 = [x for x in filenames if sample_name in x and "_R1_" in x][0]
    sample_dict[sample_name] = fastq_file_r1

# Creates necessary directories for the analysis.
mkdir("FastQC/")
mkdir("FastQC/Raw/")
mkdir("FastQC/Trim/")
mkdir("00_log/")

# The number indicates the script that has generated the file (except for 00_log).
mkdir("02_trim/")

# Run the fastqc and trimming steps.
if slow:
    # Much lower, but less memory intensive.
    sample_dict = trimming_files_slow(sample_dict, adapter, run)
else:
    eval_fastq_files(sample_dict, "FastQC/Raw", adapter, run, processes)
    sample_dict = trimming_files(sample_dict, adapter, run, processes)
    eval_fastq_files(sample_dict, "FastQC/Trim", "None", run, processes)
    get_stats_fastq_files(sample_dict, run, processes)



# If we want to operate over files stored in different folders, it is necessary to run this script on each folder separately.
# The append option allows us to append the sample dict of the new samples to the previous one.
if append:
    with open("00_log/1_2_fastq.json", "r") as jsonfile:
        appendix = sample_dict
        sample_dict = json.load(jsonfile)
        sample_dict.update(appendix)

# Writes out the sample dict as a json file in order to be exchangeable with the other scripts.
with open("00_log/1_2_fastq.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)
