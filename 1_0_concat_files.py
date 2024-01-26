import argparse
import os
import glob
import multiprocessing
from functions.libs import concatenate_files


parser = argparse.ArgumentParser()
parser.add_argument("-D", "--project")
parser.add_argument("-X", "--rem_sampls")
parser.add_argument("-R", "--run")
parser.add_argument("-P", "--processes")
args = vars(parser.parse_args())

project, rem_sampls, run, processes = (
    args["project"],
    args["rem_sampls"],
    args["run"],
    args["processes"],
)

rem_sampls = rem_sampls.split(",")

if not os.path.exists("temp_fold"):
    os.mkdir("temp_fold")

listFiles = glob.glob("{}/*_R*_001.fastq.gz".format(args["project"]))
sampleNames = list(
    set(
        [
            os.path.basename(sampleFile).split("_R1_001.fastq.gz")[0].split("_L00")[0]
            for sampleFile in listFiles
        ]
    )
)
sampleNames = [sampleName for sampleName in sampleNames if sampleName not in rem_sampls]

if processes:
    processes = int(processes)
else:
    processes = len(sampleNames)

with multiprocessing.Pool(processes) as pool:
    pool.map(
        concatenate_files,
        [(sampleName, listFiles, run) for sampleName in sampleNames],
    )

os.system("rm {}/*".format(project))
os.system("mv temp_fold/* {}".format(project))
os.system("rm -r temp_fold")
