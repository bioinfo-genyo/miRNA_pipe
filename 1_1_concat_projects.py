""" This is the second step of file preparation, to reunite all fastq files from different sequencing projects in one single folder.
Args:
    -D, --final-dir (str): Directory to store all the files.
    -X, --rem-files (str): A comma-separated list of samples to be removed.
    -P, --project (str): A comma-separated list of the projects' names.
    -L, --run (str): The run number."""

import argparse
from functions.libs import list_dir_files, mkdir, copy_files

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--final-dir")
parser.add_argument("-A", "--rem-files")
parser.add_argument("-P", "--project")
parser.add_argument("-L", "--run")
args = vars(parser.parse_args())

final_dir, rem_files, project, run = (
    args["final_dir"],
    args["rem_files"],
    args["project"],
    args["run"],
)

if final_dir != project:
    if run == "1":
        mkdir(final_dir)
        projects = project.split(",")
        for project in projects:
            lfiles = list_dir_files(project)
            [copy_files(lfile, lfile.replace(project, final_dir)) for lfile in lfiles]
