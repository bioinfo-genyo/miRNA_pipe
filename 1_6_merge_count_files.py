import argparse
from functions.libs import merge_count_files

parser = argparse.ArgumentParser()
parser.add_argument("-R", "--run")
args = vars(parser.parse_args())

run = args["run"]

# Merges all the sample counts to create the count matrix.
merge_count_files(run)
