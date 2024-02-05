import argparse
import json
from functions.libs import list_dir_files, get_sample_name

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--input-dir")
args = vars(parser.parse_args())

input_dir = args["input_dir"]

filenames = list_dir_files(input_dir, "fastq.gz")
sample_names = get_sample_name(filenames)

sample_dict = {}
for sample_name in sample_names:
    fastq_file_r1 = [x for x in filenames if sample_name in x][0]
    sample_dict[sample_name] = fastq_file_r1

with open("00_log/1_2_ref.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)
