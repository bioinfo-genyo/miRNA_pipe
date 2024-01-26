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

parser = argparse.ArgumentParser()
parser.add_argument("-R", "--ref")
parser.add_argument("-T", "--tax")
parser.add_argument("-G", "--gff")
parser.add_argument("-K", "--kegg")
parser.add_argument("-M", "--use-mirbase")
parser.add_argument("-L", "--run")
parser.add_argument("-P", "--processes")
args = vars(parser.parse_args())


with open("00_log/1_4_bam.json", "r") as jsonfile:
    sample_dict = json.load(jsonfile)

with open("00_log/1_4_mirna_counts.json", "r") as jsonfile:
    mirna_counts = json.load(jsonfile)

mkdir("04_counts")

reference_folder, tax, gff, kegg, use_mirbase, run, processes = (
    args["ref"],
    args["tax"],
    args["gff"],
    args["kegg"],
    args["use_mirbase"],
    args["run"],
    args["processes"],
)
processses = int(processes)


gtf_files = prepare_biotypes(reference_folder, gff, tax, biotypes="miRNA")
ref_file = f"{reference_folder}/mature.fa"
mirbaseDB = filter_mirbase(kegg, ref_file)

if use_mirbase != "0":
    download_file(use_mirbase, f"{reference_folder}/miRNA.gff3")
    gtf_files["miRNA"] = f"{reference_folder}/miRNA.gff3"


for biotype in gtf_files:
    gtf_file = gtf_files[biotype]
    sample_files_bio = quantify_biotype(sample_dict, gtf_file, biotype, run, processes)
    if biotype == "miRNA":
        quantify_samples(sample_files_bio, mirna_counts, run, processes)
        sample_files_bio = concat_mirna_samples(
            sample_files_bio, mirna_counts, use_mirbase, mirbaseDB, processes
        )
