""" This is the first step of file preparation, to concatenate all the fastq files for each sample is we count with lane_splitted data.
Args:
    -D, --project (str): The name of the project.
    -X, --rem_samples (str): A comma-separated list of samples to be removed.
    -R, --run (str): A run control variable.
    -P, --processes (str): The number of processes to use. If no processes are specified, use the number of samples to maximize parallelization.
"""

import argparse
import os
import glob
import multiprocessing
from functions.libs import concatenate_files


parser = argparse.ArgumentParser()
parser.add_argument("-D", "--project")
parser.add_argument("-X", "--rem_samples")
parser.add_argument("-R", "--run")
parser.add_argument("-P", "--processes", type=str, default="sample")
args = vars(parser.parse_args())

project, rem_samples, run, processes = (
    args["project"],
    args["rem_samples"],
    args["run"],
    args["processes"],
)

# Obtain a list with the samples to be removed.
rem_samples = rem_samples.split(",")

# Make a temporary folder to store the intermediate files.
if not os.path.exists("temp_fold"):
    os.mkdir("temp_fold")

# Get a list of all the fastq files in the directory.
# Filtering using the standard illumina naming.
listFiles = glob.glob("{}/*_R*_001.fastq.gz".format(args["project"]))
# Get the sample names by removing the sequencing informational suffix.
sampleNames = list(
    set(
        [
            os.path.basename(sampleFile).split("_R1_001.fastq.gz")[0].split("_L00")[0]
            for sampleFile in listFiles
        ]
    )
)
# Suppress the samples to be removed.
sampleNames = [
    sampleName for sampleName in sampleNames if sampleName not in rem_samples
]

# Concatenate the files using multiprocessing.
if processes == "sample":
    processes = len(sampleNames)

with multiprocessing.Pool(processes) as pool:
    pool.map(
        concatenate_files,
        [(sampleName, listFiles, run) for sampleName in sampleNames],
    )

# Remove temporary files and folders.
os.system("rm {}/*".format(project))
os.system("mv temp_fold/* {}".format(project))
os.system("rm -r temp_fold")
