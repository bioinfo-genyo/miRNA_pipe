""" 
This is the third step of the analysis, mapping the reads to the reference genome. We also do a initial miRNA assignment in order to filter out the reads that exactly match the miRNAs of the given database. That way those counts are directly taken into account and reduces computing time expended on the mapping step.

    Args:
        -r, --ref (str): The reference directory where all files will be stored.
        -K, --kegg (str): The kegg biotype to filter de miRNA database (hsa for human).
        -R, --run (str): Run control variable (True to run).
        -T, --threads (int): The number of threads to use in applications that allow multithreading. Default is the number of CPU threads.
        -P, --processes (str): The number of cpu threads to use. Default is 4. If 0 is specified, use the number of samples to maximize parallelization.
"""

import argparse
import json
from functions.libs import (
    mkdir,
    download_file,
    filter_mirbase,
    mirbase_sequence_assign,
    align_samples,
    quality_mapping_samples,
    create_sample_dict,
)
from multiprocessing import cpu_count

# Gets the command line arguments with argparse.
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", type=str)
parser.add_argument("-K", "--kegg", type=str)
parser.add_argument("-T", "--threads", type=int, default=cpu_count())
parser.add_argument("-P", "--processes", type=int, default=4)
parser.add_argument("-R", "--run", type=bool, default=False)
args = vars(parser.parse_args())

# Assign the command line arguments to variables.
reference_folder, kegg, threads, processes, run = (
    args["ref"],
    args["kegg"],
    args["threads"],
    args["processes"],
    args["run"],
)

# Gets reference location.
with open("00_log/1_3_ref.json", "r") as jsonfile:
    bowtie_reference = json.load(jsonfile)

try:
    # Loads the sample dictionary.
    with open("00_log/1_2_fastq.json", "r") as jsonfile:
        sample_dict = json.load(jsonfile)
except FileNotFoundError:
    # If not file is found, builds the sample dictionary from the 1_2_fastqc_trimming.py output files.
    sample_dict = create_sample_dict("02_trim/", "_trimmed.fastq.gz")

# The number indicates the script that has generated the file (except for 00_log).
mkdir("04_bam/")

# Downloads and filters the mirBase database (keep only our biotype of interest, human for example)
download_file("https://mirbase.org/download/mature.fa", f"{reference_folder}/mature.fa")
ref_file = f"{reference_folder}/mature.fa"
mirbaseDB = filter_mirbase(kegg, ref_file)

# Performs read pre-counting.
sample_dict, mirna_counts = mirbase_sequence_assign(sample_dict, mirbaseDB, processes)

# Aligns not counted reads.
sample_dict = align_samples(sample_dict, bowtie_reference, threads, processes, run)
quality_mapping_samples(sample_dict, mirna_counts, processes, run)

# Writes out the sample dict as a json file in order to be exchangeable with the other scripts.
with open("00_log/1_4_bam.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)

# Writes out the path to the counts as a json file in order to be exchangeable with the other scripts.
with open("00_log/1_4_mirna_counts.json", "w") as jsonfile:
    json.dump(mirna_counts, jsonfile, indent=4)
