import argparse
import json
from functions.libs import list_dir_files,get_sample_name,mkdir,eval_fastq_files,trimming_files,get_stats_fastq_files

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--input-dir")
parser.add_argument("-A", "--adapter")
parser.add_argument("-L", "--run")
args = vars(parser.parse_args())

input_dir = args["input_dir"]
adapter = args["adapter"]
run = args["run"]

filenames = list_dir_files(input_dir,"fastq.gz")
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

eval_fastq_files(sample_dict,"FastQC/Raw",adapter,run)
sample_dict = trimming_files(sample_dict,adapter,run)
eval_fastq_files(sample_dict,"FastQC/Trim","None",run)
get_stats_fastq_files(sample_dict,run)

with open("00_log/1_2_fastq.json","w") as jsonfile:
    json.dump(sample_dict,jsonfile,indent=4)