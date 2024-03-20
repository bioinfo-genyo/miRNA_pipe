import argparse
import json
import os
from functions.libs import download_file, prepare_ref

parser = argparse.ArgumentParser()
parser.add_argument("-F", "--fasta")
parser.add_argument("-R", "--ref")
args = vars(parser.parse_args())

fasta, ref = (args["fasta"], args["ref"])

destfile = f"{ref}/{os.path.basename(fasta)}"
download_file(fasta, destfile)
prepare_ref(destfile, ref)

reference = {"index": f"{ref}/Bowtie/genome", "fasta": destfile}
with open("00_log/1_3_ref.json", "w") as jsonfile:
    json.dump(reference, jsonfile, indent=4)
