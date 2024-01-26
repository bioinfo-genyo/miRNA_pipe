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

with open("00_log/1_3_ref.json", "r") as jsonfile:
    reference = json.load(jsonfile)

with open("00_log/1_2_fastq.json", "r") as jsonfile:
    sample_dict = json.load(jsonfile)

mkdir("03_bam")

reference_folder, kegg, run, processes = (
    args["ref"],
    args["kegg"],
    args["run"],
    args["processes"],
)
processes = int(processes)

download_file("https://mirbase.org/download/mature.fa", f"{reference_folder}/mature.fa")
ref_file = f"{reference_folder}/mature.fa"
mirbaseDB = filter_mirbase(kegg, ref_file)

sample_dict, mirna_counts = mirbase_sequence_assign(sample_dict, mirbaseDB, processes)
sample_dict = align_samples(sample_dict, reference, run, processes)
quality_mapping_samples(sample_dict, mirna_counts, run, processes)

with open("00_log/1_4_bam.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)

with open("00_log/1_4_mirna_counts.json", "w") as jsonfile:
    json.dump(mirna_counts, jsonfile, indent=4)
