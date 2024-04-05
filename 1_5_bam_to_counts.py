"""
This script is the fifth step of the analysis. We will count the aligned reads to the miRNAs of the given database using featureCounts. After that, the pre-counted reads will be added to the results of featureCounts, obtaining the final count files for each sample.

    Args:
        -R, --ref (str): The genome reference directory where all files will be stored.
        -T, --tax (str): The taxonomic code to use (9606 for human).
        -G, --gff (str): The url to the file of the genome annotation gff file (input for featureCounts). 
        -K, --kegg (str): The kegg biotype to filter de miRNA database (hsa for human).
        -M, --use-mirbase (str): URL to the mirbase database to download it. In this case is also used as a flag. When the species we are analyzing is not in mirBase, left use_mirbase emtpy and the data from the available database will be used
        -P, --processes (str): The number of cpu threads to use. Default is 4. If 0 is specified, use the number of samples to maximize parallelization.
"""

import argparse
import json
from functions.libs import (
    mkdir,
    prepare_biotypes,
    filter_mirbase,
    download_file,
    concat_mirna_samples,
    quantify_biotype,
    quantify_samples,
)

# Gets the command line arguments with argparse.
parser = argparse.ArgumentParser()
parser.add_argument("-R", "--ref", type=str)
parser.add_argument("-T", "--tax", type=str)
parser.add_argument("-G", "--gff", type=str)
parser.add_argument("-K", "--kegg", type=str)
parser.add_argument("-M", "--use-mirbase", type=str, default=None)
parser.add_argument("-L", "--run", type=bool, default=False)
parser.add_argument("-P", "--processes", type=int, default=4)
args = vars(parser.parse_args())

# Assign the command line arguments to variables.
reference_folder, tax, gff, kegg, use_mirbase, run, processes = (
    args["ref"],
    args["tax"],
    args["gff"],
    args["kegg"],
    args["use_mirbase"],
    args["run"],
    args["processes"],
)

# Loads the sample dictionary and the paths to count files.
with open("00_log/1_4_bam.json", "r") as jsonfile:
    sample_dict = json.load(jsonfile)

with open("00_log/1_4_mirna_counts.json", "r") as jsonfile:
    mirna_counts = json.load(jsonfile)

# The number indicates the script that has generated the file (except for 00_log).
mkdir("05_counts/")

# Filters the mirBase database (keep only our species of interest, human for example).
gtf_files = prepare_biotypes(reference_folder, gff, tax, biotypes="miRNA")
ref_file = f"{reference_folder}/mature.fa"
mirbaseDB = filter_mirbase(kegg, ref_file)

if use_mirbase != "0":
    download_file(use_mirbase, f"{reference_folder}/miRNA.gff3")
    gtf_files["miRNA"] = f"{reference_folder}/miRNA.gff3"

# Quantifies only the miRNAs.
for biotype in gtf_files:
    gtf_file = gtf_files[biotype]
    sample_files_bio = quantify_biotype(sample_dict, gtf_file, biotype, run, processes)
    if biotype == "miRNA":
        quantify_samples(sample_files_bio, mirna_counts, run, processes)
        sample_files_bio = concat_mirna_samples(
            sample_files_bio, mirna_counts, use_mirbase, mirbaseDB, processes
        )
