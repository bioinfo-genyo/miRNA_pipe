""" 
In this script we prepare the reference genome for the alignment.

    Args:
        -F, --fasta (str): The url to the file of the reference genome.
        -R, --ref (str): The reference directory where all file will be stored. Default is current directory.
"""

import argparse
import json
import os
from functions.libs import download_file, prepare_ref

# Gets the command line arguments with argparse.
parser = argparse.ArgumentParser()
parser.add_argument("-F", "--fasta", type=str)
parser.add_argument("-R", "--ref", type=str, default=".")

# Assign the command line arguments to variables.
args = vars(parser.parse_args())
fasta, ref = (args["fasta"], args["ref"])

# Downloads the file into the reference folder. The basename of the url will be the name of the fasta file in the file system.
destfile = f"{ref}/{os.path.basename(fasta)}"
download_file(fasta, destfile)
prepare_ref(destfile, ref)

# Writes out the reference location as a json file in order to be exchangeable with the other scripts.
# WARNING: /Bowtie/genome indicates a basename for all de reference files, not a directory, do not add / at the end.
reference = {"index": f"{ref}/Bowtie/genome", "fasta": destfile}
with open("00_log/1_3_ref.json", "w") as jsonfile:
    json.dump(reference, jsonfile, indent=4)
