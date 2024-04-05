""" 
This is the third step of the analysis, mapping the reads to the reference genome. We also do a initial miRNA assignment in order to filter out the reads that exactly match the miRNAs of the given database. That way those counts are directly taken into account and reduces computing time expend on the mapping step.

    Args:
        -R, --ref (str): The reference directory where all files will be stored.
        -K, --kegg (str): The kegg biotype to filter de miRNA database (hsa for human).
        -L, --run (str): Run control variable (1 to run).
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
)

# Gets the command line arguments with argparse.
parser = argparse.ArgumentParser()
parser.add_argument("-R", "--ref", type=str)
parser.add_argument("-K", "--kegg", type=str)
parser.add_argument("-L", "--run", type=bool, default=False)
parser.add_argument("-P", "--processes", type=int, default=4)
args = vars(parser.parse_args())

# Assign the command line arguments to variables.
reference_folder, kegg, run, processes = (
    args["ref"],
    args["kegg"],
    args["run"],
    args["processes"],
)

# Gets reference location.
with open("00_log/1_3_ref.json", "r") as jsonfile:
    bowtie_reference = json.load(jsonfile)

# Loads sample dict.
with open("00_log/1_2_fastq.json", "r") as jsonfile:
    sample_dict = json.load(jsonfile)

# The number indicates the script that has generated the file (except for 00_log).
mkdir("04_bam/")

# Downloads and filters the mirBase database (keep only our biotype of interest, human for example)
download_file("https://mirbase.org/download/mature.fa", f"{reference_folder}/mature.fa")
ref_file = f"{reference_folder}/mature.fa"
mirbaseDB = filter_mirbase(kegg, ref_file)

# Performs read pre-counting.
sample_dict, mirna_counts = mirbase_sequence_assign(sample_dict, mirbaseDB, processes)

# Aligns not counted reads.
sample_dict = align_samples(sample_dict, bowtie_reference, run, processes)
quality_mapping_samples(sample_dict, mirna_counts, run, processes)

# Writes out the sample dict as a json file in order to be exchangeable with the other scripts.
with open("00_log/1_4_bam.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)

# Writes out the path to the counts as a json file in order to be exchangeable with the other scripts.
with open("00_log/1_4_mirna_counts.json", "w") as jsonfile:
    json.dump(mirna_counts, jsonfile, indent=4)
