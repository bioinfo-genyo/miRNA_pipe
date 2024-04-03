""" 
This is the third step of the analysis, mapping the reads to the reference genome. We also do a initial miRNA assignment in order to filter out the reads that exactly match the miRNAs of the given database. That way those counts are directly taken into account and reduces computing time expend on the mapping step.

    Args:
        -R, --ref (str): The reference directory where all files will be stored.
        -K, --kegg (str): The url to the file of the kegg database.
        -L, --run (str): Run control variable (1 to run).
        -P, --processes (str): The number of cpu threads to use. If no number of threads is specified, use the number of samples to maximize parallelization.
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

parser = argparse.ArgumentParser()
parser.add_argument("-R", "--ref")
parser.add_argument("-K", "--kegg")
parser.add_argument("-L", "--run")
parser.add_argument("-P", "--processes")
args = vars(parser.parse_args())

# Gets reference location.
with open("00_log/1_3_ref.json", "r") as jsonfile:
    bowtie_reference = json.load(jsonfile)

# Loads sample dict.
with open("00_log/1_2_fastq.json", "r") as jsonfile:
    sample_dict = json.load(jsonfile)

mkdir("03_bam")

reference_folder, kegg, run, processes = (
    args["ref"],
    args["kegg"],
    args["run"],
    args["processes"],
)

if processes:
    processes = int(processes)

# Downloads and filters the mirBase database (keep only our biotype of interest, human for example)
download_file("https://mirbase.org/download/mature.fa", f"{reference_folder}/mature.fa")
ref_file = f"{reference_folder}/mature.fa"
mirbaseDB = filter_mirbase(kegg, ref_file)

# Performs read pre-counting.
sample_dict, mirna_counts = mirbase_sequence_assign(sample_dict, mirbaseDB, processes)

# Aligns not counted reads.
sample_dict = align_samples(sample_dict, bowtie_reference, run, processes)
quality_mapping_samples(sample_dict, mirna_counts, run, processes)

with open("00_log/1_4_bam.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)

with open("00_log/1_4_mirna_counts.json", "w") as jsonfile:
    json.dump(mirna_counts, jsonfile, indent=4)
