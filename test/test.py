import json
from random import sample

sample_dict = {"sample1": "fastq1.fastq", "sample2": "fastq2.fastq"}

with open("1_2_fastq.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)

sample_dict = {"sample3": "fastq3.fastq", "sample4": "fastq4.fastq"}

with open("1_2_fastq.json", "r") as jsonfile:
    appendix = sample_dict
    sample_dict = json.load(jsonfile)
    sample_dict.update(appendix)

with open("1_2_fastq.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)
