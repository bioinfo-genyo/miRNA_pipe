#### IMPORT LIBRARIES ####

from typing import Generator, Literal
from multiprocessing import Pool, cpu_count
from subprocess import run as bash  # Avoids namespace conflicts.
import collections
import shutil
import os
import glob
import tqdm
import requests
import gzip
import numpy as np
import pandas as pd
import re


#### SET CPU THREADS ####

num_threads = cpu_count()


#### FUNCTIONS ####


def detect_paired_single(
    sampleName: str, listFiles: list
) -> Literal["paired", "single"]:
    """
    Detects if the sample is paired or single-end based on the given sample name and list of files.
    Args:
        sampleName (str): The name of the sample to detect.
        listFiles (list): The list of files to search for the sample.
    Returns:
        str: "paired" if the sample has 2 matching files, "single" otherwise.
    """

    sampleFiles = [
        sampleFile
        for sampleFile in listFiles
        if sampleName in sampleFile
        and "_L001_" in sampleFile  # Detects illumina standard file naming.
    ]
    if len(sampleFiles) == 2:
        return "paired"
    else:
        return "single"


def shutil_python(output_file: str, input_files: list[str]) -> None:
    """
    Concatenates the data from multiple input files and writes the result to the specified output file.
    Args:
        output_file (str): The path to the output file.
        input_files (List[str]): The list of paths to the input files.
    Returns:
        None
    """

    # Concatenate the data from the input files
    with open(output_file, "wb") as out:
        for input_file in input_files:
            with open(input_file, "rb") as f_in:
                shutil.copyfileobj(f_in, out)
    # Compress the concatenated data
    bash(["gzip", output_file], check=True)


def zcat_files(output_file: str, input_files: list[str]) -> None:
    """
    Compresses and concatenates the input files into the specified output file.

    Args:
        output_file (str): The name of the output file to be created.
        input_files (list[str]): The list of input files to be concatenated and compressed.

    Returns:
        None
    """

    os.system("zcat {} >{}".format(" ".join(input_files), output_file))
    os.system("gzip {}".format(output_file))


def concatenate_files(
    sampleName: str, listFiles: list[str], run: bool = False
) -> dict[any, dict[any, any]]:
    """
    Concatenates files based on sample name and list of files and returns a dictionary
    containing sample name and corresponding files.
    This function is thought to concatenate lane-split files.

    Args:
        sampleName (str): The name of the sample to concatenate.
        listFiles (list): The list of files to concatenate.
        run (bool, optional): Whether to run the concatenation. Defaults to False.

    Returns:
        dict: A dictionary with sample name as key and a nested dictionary containing
        file information as value.
    """

    paired_or_single = detect_paired_single(sampleName, listFiles)
    # Create the strands list to add to the sample name.
    if paired_or_single == "paired":
        strands = ["R1", "R2"]
    else:
        strands = ["R1"]
    sampleDict = {}
    # Change the sample name to the required format (remove "-" and replace with "_" for R compatibility)
    sampleName_change = sampleName.replace("-", "_")
    # List all the files associated with the sample name.
    sampleFiles = [sampleFile for sampleFile in listFiles if sampleName in sampleFile]
    # Concatenate the files of different strands separated.
    for strand in strands:
        input_files = [sampleFile for sampleFile in sampleFiles if strand in sampleFile]
        input_files.sort()
        output_file = "temp_fold/" + os.path.basename(
            input_files[0].replace(
                "_L001_", "_"
            )  # Names the output file as the input files without the lane information.
        )
        output_file = output_file.replace("-", "_").replace(
            ".gz", ""
        )  # Replaces the "-" and ".gz" with "_" in the name.
        if run:
            if os.path.exists(output_file + ".gz"):
                os.remove(output_file + ".gz")
            zcat_files(output_file, input_files)  # Concatenate and compress the files.
        sampleDict[strand] = output_file + ".gz"
    return {sampleName_change: sampleDict}


def mkdir(dir: str) -> None:
    """
    Function to create a directory
    """

    if not os.path.exists(dir):
        os.mkdir(dir)


def list_dir_files(dir: str, pattern: str = None) -> list[str]:
    """
    Function to list the files of a directory. If pattern is specified, only files matching the pattern will be returned.
    """

    if not pattern:
        files = glob.glob(f"{dir}/*")
    else:
        files = glob.glob(f"{dir}/*{pattern}*")
    return files


def copy_files(src_file: str, dst_file: str) -> None:
    """
    Function to copy files from source to dest
    """

    shutil.copy(src_file, dst_file)


def rm_file(file: str) -> None:
    """
    Function to remove file if it exists
    """

    if os.path.exists(file):
        os.remove(file)


def download_file(url: str, filename: str, force: bool = False) -> None:
    """
    Function to download data files. If the file already exists, it will not be downloaded again.

    Args:
        url (str): The URL of the file to download.
        filename (str): The name of the file to download.
        force (bool, optional): Whether to force download even if file exists. Defaults to False.

    Returns:
        None
    """

    # First check if the file already exists.
    if (
        not os.path.exists(filename)
        and not os.path.exists(filename.replace(".gz", ""))
        and not os.path.exists(filename + ".gz")
    ) or force:  # Force download.
        rm_file(filename)
        rm_file(filename.replace(".gz", ""))
        with open(filename, "wb") as f:
            # request is a library for making HTTP requests.
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                # Get progress information for tqdm.
                total = int(r.headers.get("content-length", 0))
                # tqdm is a progress verbose (progress bar, performing step...) library.
                # tqdm has many interesting parameters. Feel free to experiment!
                tqdm_params = {
                    "desc": url,
                    "total": total,
                    "miniters": 1,
                    "unit": "B",
                    "unit_scale": True,
                    "unit_divisor": 1024,
                }
                with tqdm.tqdm(**tqdm_params) as pb:
                    for chunk in r.iter_content(chunk_size=8192):
                        pb.update(len(chunk))
                        f.write(chunk)


def get_sample_name(file_names: list[str]) -> list[str]:
    """
    Function to list the sample names of a list of fastq files.
    """

    return list(
        set(
            [
                os.path.basename(file).split("_R1_")[0]
                for file in file_names
                if "_R1_" in file
            ]
        )
    )


def read_gzfile(filename: str) -> Generator[bytes | str, any, None]:
    # The use of a generator allows better memory management (loads the file in chunks).
    """
    Function to read a .gz file
    """

    with gzip.open(filename, "rt") as f:
        for line in f:
            yield line.rstrip()


def write_log(logfile: str, text: str, mode: str) -> None:
    """
    Writes the given text to the specified log file using the provided mode.

    Args:
        logfile (str): The path to the log file.
        text (str): The text to be written to the log file.
        mode (str): The mode in which the file should be opened.

    Returns:
        None
    """

    with open(logfile, mode) as write_file:
        write_file.write(text)


def eval_fastq_file(args: tuple) -> None:
    """
    A function to evaluate a fastq file. Performs fastq and gets the number of reads from the file.

    Args:
        args (tuple): The arguments must be input in a single tuple with the following components (allows multiprocessing):
            sample_path (str): The path to the fastq file.
            output (str): The output directory for the results.
            adapter (str): The adapter to be used.
            threads (int): The number of threads to be used.
            run (bool): Whether to run the evaluation.

    Returns:
        None
    """

    # Unpack the arguments. The first argument is the sample path, the second is the output directory, the third is the adapter, the fourth is the number of threads, and the fifth is the run flag. This allows to accept a tuple as input, making the function compatible with multiprocessing.
    sample_name, sample_path, output, adapter, threads, run = args

    if run:
        # Run fastqc.
        bash(f"fastqc {sample_path} -o {output} -t {threads}", shell=True)
        # Get the filename (! filename != sample_name).
        filename = os.path.basename(sample_path)
        # Create the path for the log file.
        log_file = f"00_log/{sample_name}.log"

        # If  the adapter is included, that means the fastq file is raw.
        if adapter:
            mode = "w"
            text = "############ RAW READS ##############\n\n"
            write_log(log_file, text, mode)
            mode = "a"
            # Creates a list of lines from the fastq file. The file is contained in a generator that will only read the lines when requested.
            lines = read_gzfile(sample_path)
            # Because the fastq file has 4 lines per sequence, the numpy array is used to split the list into chunks of 4 lines with the zip function.
            gzip_cont = np.array(list(zip(*[lines] * 4)))
            # Write the number of reads to the log file.
            text = f"{filename} has {int(len(gzip_cont))} reads\n\n"
            write_log(log_file, text, mode)

        # If no adapter is specified, the fastq file is trimmed.
        else:
            mode = "a"
            text = "############ TRIM READS ##############\n\n"
            write_log(log_file, text, mode)
            mode = "a"
            # Creates a list of lines from the fastq file. The file is contained in a generator that will only read the lines when requested.
            lines = read_gzfile(sample_path)
            # Because the fastq file has 4 lines per sequence, the numpy array is used to split the list into chunks of 4 lines with the zip function.
            gzip_cont = np.array(list(zip(*[lines] * 4)))
            # Write the number of reads to the log file.
            text = f"{filename} has {int(len(gzip_cont))} reads\n\n"
            write_log(log_file, text, mode)


def eval_fastq_files(
    sample_dict: dict,
    output: str,
    adapter: str = None,
    threads: int = num_threads,
    processes: int = 4,
    run: bool = False,
) -> None:
    """
    Runs eval_fastq_file in parallel using multiprocessing and performs multiqc.

    Args:
        sample_dict (dict): A dictionary containing sample information.
        output (str): The output directory for the evaluated files.
        adapter (str): The adapter to be used for evaluation.
        threads (int, optional): The number of threads to be used. Defaults to total number of CPU threads.
        run (bool, optional): Whether to run the evaluation. Defaults to False.
        processes (str, optional): The number of processes to be used for multiprocessing. Defaults to the number of samples to maximize parallelization.

    Returns:
        None: This function does not return anything.
    """

    if processes == 0:
        processes = len(sample_dict)

    with Pool(processes) as pool:
        pool.map(
            eval_fastq_file,
            [
                # Creates a tuple for each sample to be used as input for the eval_fastq_file function.
                (sample_name, sample_dict[sample_name], output, adapter, threads, run)
                for sample_name in sample_dict
            ],
        )

    # Run multiqc.
    bash(f"multiqc {output} --outdir {output} --force", shell=True)


def remove_umi_delete_adapter(fastq_file: str, adapter: str, outfile: str) -> int:
    """
    Remove UMIs and delete adapters from a given fastq file and save the processed data to an output file.

    Args:
        fastq_file (str): The input fastq file path.
        adapter (str): The adapter sequence to be removed. WARNING: If the adapter sequence is not correct, the function will output an empty file.
        outfile (str): The output file path to save the processed data.

    Returns:
        int: The total count of duplicated lines after processing.
    """

    # Creates a list of lines from the fastq file. The file is contained in a generator that will only read the lines when requested.
    lines_list = read_gzfile(fastq_file)
    # Because the fastq file has 4 lines per sequence, the numpy array is used to split the list into chunks of 4 lines with the zip function.
    gzip_cont = np.array(list(zip(*[lines_list] * 4)))
    # Create a variable to store the unique sequences.
    unique_elements = set()
    # Create a variable to store the processed reads.
    filtered_lines = []
    # Create a variable to store the duplicated reads for logging purposes.
    duplicated_lines = []
    # Loop through each sequence and check if the adapter is in the sequence. If it is, remove the adapter (+ the UMI) from the sequence and add the new sequence and quality to the list of processed reads.
    for sublist in gzip_cont:
        # Unpack the sequence and quality from the sublist.
        seq, quality = sublist[1], sublist[3]
        if seq not in unique_elements:
            # WARNING: If the adapter sequence is not correct, the function will output an empty file.
            if adapter in seq:
                unique_elements.add(seq)
                # Find the adapter start position in the read.
                idx_pos = seq.rfind(adapter)
                # Remote all bases from that position to the end.
                seq_wo_adapter = seq[:idx_pos]
                quality_wo_adapter = quality[:idx_pos]
                # Add the new sequence and quality to the list of processed reads.
                sublist[1] = seq_wo_adapter
                sublist[3] = quality_wo_adapter
                filtered_lines.append(sublist)
        else:
            duplicated_lines.append(seq)
    # Counts repeated reads. The unique_dup list contains the unique set of duplicated reads to take into account the ones that have been processed.
    unique_dup = list(set(duplicated_lines))
    duplicated_lines = len(duplicated_lines) + len(unique_dup)
    # Creates a new gzip file with the processed un-duplicated reads.
    gzip_cont = "\n".join(["\n".join(seq) for seq in filtered_lines])
    with gzip.open(outfile, "wb") as f:
        f.write(gzip_cont.encode())
    return duplicated_lines


def run_trimming(args) -> dict[str, str]:
    """
    Trims fastq files using the specified adapter and run settings.

    Args:
        args (tuple): The arguments must be input in a single tuple with the following components (allows multiprocessing):
            sample_name (str): The sample name.
            fastq_file (str): The fastq file path.
            adapter (str): The adapter sequence for trimming.
            threads (int): The number of threads to use for trimming.
            run (bool): Whether to run the trimming.

    Returns:
        Sample dictionary with sample name as key and fastq file path as value.
    """

    # Unpack the arguments. The first argument is the sample name, the second is the fastq file path, the third is the adapter sequence, the fourth is the number of threads to use, and the fifth is the run flag. This allows to accept a tuple as input, making the function compatible with multiprocessing.
    sample_name, fastq_file, adapter, threads, run = args

    # For each step, the file name is updated to indicate the process performed over the sample.
    # This file name change is registered in a variable called sample_dict. This is a dictionary with sample_name as key and fastq_file path as a value.
    # That variable is the key for the integration of all the functions in the pipeline, as it allows the correct localization of the appropriate file for each step.
    outFileUMI = f"02_trim/{sample_name}_umi.fastq.gz"
    outFileCut = f"02_trim/{sample_name}_trimmed.fastq.gz"
    if run:
        print("Running umi removal {}".format(fastq_file))
        duplicated_lines = remove_umi_delete_adapter(fastq_file, adapter, outFileUMI)

        # Updates the log for that sample.
        mode = "a"
        text = "############ DUPS READS ##############\n\n"
        log_file = f"00_log/{sample_name}.log"
        write_log(log_file, text, mode)
        mode = "a"
        text = f"{sample_name} has {duplicated_lines} duplicated reads\n\n"
        write_log(log_file, text, mode)

        # Cutadapt is used for quality trimming
        print("Running cutadapt {}".format(fastq_file))
        bash(
            f"cutadapt --quiet -j {threads} -m 10 -M 40 -q 10 {outFileUMI} -o {outFileCut}",
            shell=True,
        )

    # Remove the intermediate files
    rm_file(outFileUMI)
    # Returns a dict containing the sample name as key and the trimmed fastq file as value.
    # This would be the basic return value for the functions in order to build the sample dict.
    # This dict relates the sample names (the identifier the functions work with) with the path of the file in order to be find by the function calls.
    # After each step, this dict is updated with the new file names for each samples as processing steps are performed over them.
    return {sample_name: outFileCut}


def trimming_files(
    sample_dict: dict,
    adapter: str,
    threads: int = num_threads,
    processes: int = 4,
    run: bool = False,
) -> dict[str, str]:
    """
    Runs run_trimming in parallel using multiprocessing and returns a dictionary containing the trimmed files.

    Args:
        sample_dict (dict): A dictionary containing the samples to be trimmed.
        adapter (str): The adapter used for trimming.
        run (int): The run number.
        threads (int): The number of threads to use for trimming. Defaults to total number of CPU threads.
        processes (str, optional): The number of processes to use for multiprocessing. Defaults to "sample".
        run (bool, optional): Whether to run the trimming. Defaults to False.

    Returns:
        dict: A dictionary containing the trimmed files.
    """

    if processes == 0:
        processes = len(sample_dict)

    with Pool(processes) as pool:
        sample_dict = pool.map(
            run_trimming,
            [
                # Creates a tuple for each sample to be used as input for the run_trimming function.
                (sample_name, sample_dict[sample_name], adapter, threads, run)
                for sample_name in sample_dict
            ],
        )

    # Collects all the generated sample_dict outputs from each call to run_trimming and creates a single dictionary where the keys are the sample name and the values the file paths for each sample. For each step, the file name is updated to indicate the performed process. In this case, teh suffix "_trimmed is added".
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict


def trimming_files_slow(
    sample_dict: str, adapter: str, threads: int = num_threads, run: bool = False
) -> dict[str, str]:
    """
    Trims files in the sample_dict using the specified adapter and run settings without using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing sample names as keys and sample files paths as values.
        adapter (str): The adapter used for trimming.
        threads (int, optional): The number of threads to use for trimming. Defaults to total number of CPU threads.
        run (bool, optional): Whether to run the trimming. Defaults to False.

    Returns:
        dict: A dictionary containing trimmed sample names as keys and trimmed file paths as values.
    """

    # Im case the RAM is not enough, only one sample is processed at a time.
    trimmed_dict = {}
    for sample_name in sample_dict:
        trimmed_dict.update(
            run_trimming(
                args=(sample_name, sample_dict[sample_name], adapter, threads, run)
            )
        )

    return trimmed_dict


def convert_quality_to_numeric(quality_str: str) -> list[int]:
    """
    Convert ASCII quality scores to numeric values
    """

    # Convert ASCII quality scores to numeric values
    # ord obtains the ASCII value of a character.
    quality_str_num = [ord(str(char)) - 33 for char in quality_str]
    return quality_str_num


def get_fastq_stats(args: tuple) -> None:
    """
    Calculate statistics for fastq files.

    Args:
        args (tuple): The arguments must be input in a single tuple with the following components (allows multiprocessing):
            sample_name (str): The sample name.
            fastq (str): The path to the fastq file.
            run (bool): Whether to run the evaluation.

    Returns:
        None
    """

    # Unpacks the tuple. The first element is the sample name, the second element is the fastq file, and the third element is the run flag. This allows to accept a tuple as input, making the function compatible with multiprocessing.
    sample_name, fastq, run = args

    if run:

        # Follows the same strategy for file reading as in eval_fastq_file.
        lines = read_gzfile(fastq)
        gzip_cont = np.array(list(zip(*[lines] * 4)))
        quali_seq = [line[3].rstrip() for line in gzip_cont]
        read_quality = [convert_quality_to_numeric(seq) for seq in quali_seq]

        # Obtains the maximum length of a read.
        read_lengths = [len(qual) for qual in read_quality]
        max_i = max(read_lengths)

        # Creates a dictionary to store teh quality and length statistics for each sample.
        quality_stats = {}
        lengths_stats = {}

        # Calculates the mean, median, and standard deviation of the quality scores for each position.
        # It uses numpy arrays to perform the calculations quickly.
        for pos in range(1, max_i + 1):  # +1 to include the last position.
            qual = np.array([seq[pos - 1] for seq in read_quality if len(seq) >= pos])
            mean_qual = np.mean(qual)
            median_qual = np.median(qual)
            sd_qual = np.std(qual)
            quality_stats[str(pos)] = {
                "nreads": len(qual),  # number of reads.
                "median": median_qual,  # median of the quality scores.
                "mean": mean_qual,  # mean of the quality scores.
                "sd": sd_qual,  # standard deviation of the quality scores.
            }
            lengths_stats[str(pos)] = read_lengths.count(pos)

        # Write the statistics to the log file.
        logfile = f"00_log/{sample_name}.log"
        mode = "a"
        ### quality stats
        text = "############ QUALITY READS ##############\n\n"
        write_log(logfile, text, mode)
        text = "bp\tmean\n"
        write_log(logfile, text, mode)
        for i in quality_stats:
            text = f"{i}\t{quality_stats[i]['mean']}\n"
            write_log(logfile, text, mode)

        ### lengths stats
        text = "\n############ LEN READS ##############\n\n"
        write_log(logfile, text, mode)
        text = "len\tdensity\n"
        write_log(logfile, text, mode)
        for i in lengths_stats:
            text = f"{i}\t{lengths_stats[i]}\n"
            write_log(logfile, text, mode)


def get_stats_fastq_files(
    sample_dict: dict,
    processes: int = 4,
    run: bool = False,
) -> None:
    """
    Runs get_fastq_stats in parallel using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing sample names as keys and fastq file paths as values.
        processes (int, optional): The number of processes to be used. Defaults to 4.
        run (bool, optional): Whether to run the evaluation. Defaults to False.

    Returns:
        None
    """

    if processes == 0:
        processes = len(sample_dict)

    with Pool(processes) as pool:
        pool.map(
            get_fastq_stats,
            [
                # Creates a tuple for each sample to be used as input for the get_stats_fastq_files function.
                (sample_name, sample_dict[sample_name], run)
                for sample_name in sample_dict
            ],
        )


def prepare_ref(fasta: str, ref: str, threads: int = num_threads) -> None:
    """
    Prepares reference genome index for mapping using bowtie-build.

    Args:
        fasta (str): The path to the fasta file.
        ref (str): The reference directory.

    Returns:
        None
    """

    mkdir(f"{ref}/Bowtie")
    bw_files = list_dir_files(f"{ref}/Bowtie")
    if len(bw_files) == 0:
        ########### Launch Reference ###############
        bash(f"gunzip -kf {fasta}", shell=True)
        bash(
            f"bowtie-build {fasta.replace('.gz','')} {ref}/Bowtie/genome --threads {threads}",
            shell=True,
        )
        bash(f"rm {fasta.replace('.gz','')}", shell=True)


def filter_gff(
    gene_loc: list, biotype: str, save_path: str, header: str, idmapcont_ndict: dict
) -> str:
    """
    Filters gff file by biotype and save the filtered content to a new file.

    Args:
        gene_loc (list): List of lines from the gff file.
        biotype (str): The biotype to filter the gff file by.
        save_path (str): The path to save the filtered content.
        header (str): The header to prepend to the filtered content.
        idmapcont_ndict (dict): Dictionary containing the mapping of IDs to new names.

    Returns:
        str: The path to the saved filtered content file.
    """

    # Filter biotype
    gene_loc_biotype = [
        line
        for line in gene_loc
        if line.split("\t")[-1].split("type=")[1].split(";")[0] == biotype
    ]
    # List to dictionary. miRNA name as key.
    gene_loc_biotype_dict = {
        line.split("\t")[-1].split("Name=")[1].split("_")[0]: line
        for line in gene_loc_biotype
    }
    # Filter ID

    # Get only the IDs from miRNAs present in the gft file.
    idmapcont_ndict = idmapcont_ndict[biotype]
    idmapcont_ndict = {
        key: idmapcont_ndict[key]
        for key in idmapcont_ndict
        if key in gene_loc_biotype_dict
    }

    # Only miRNAs that appear in RNAcentral.
    gene_loc_biotype_dict = {
        key: gene_loc_biotype_dict[key]
        for key in gene_loc_biotype_dict
        if key in idmapcont_ndict
    }
    # Remove duplicates and format the gff file.
    gene_loc_biotype = []
    for key in gene_loc_biotype_dict:
        line = gene_loc_biotype_dict[key].split("\t")
        line[2] = biotype  # Correct biotype.
        line[-1] = "Name={}".format(key)  # Clean last line. Only the name.
        line = "\t".join(line)
        if line not in gene_loc_biotype:
            gene_loc_biotype.append(line)

    dict_map = []
    for line in gene_loc_biotype:
        line = line.split("\t")
        name = line[-1].replace("Name=", "")
        new_name = idmapcont_ndict[name]
        line[-1] = "Name={}".format(new_name)
        line = "\t".join(line)
        dict_map.append(line)

    gtf_cont = "\n".join(dict_map)
    gtf_cont = f"{header}{gtf_cont}"
    outfile = save_path.replace(".gff3.gz", "_{}.gff3".format(biotype))
    with open(outfile, "w") as out:
        out.write(gtf_cont)
    return outfile


def prepare_biotypes(
    reference_folder: str, gff: str, tax: str, biotypes: str = "miRNA"
) -> dict[str, str]:
    """
    Prepares biotypes from a given GFF file for a specific taxonomy and return the filtered GFF files.

    Args:
        reference_folder (str): The path to the reference folder.
        gff (str): The path to the GFF file.
        tax (str): The taxonomy ID.
        biotypes (str, optional): The biotypes to prepare. Defaults to "miRNA".

    Returns:
        dict: A dictionary containing filtered GFF files for the specified biotypes.
    """

    # Downloads the gff file and saves it in the reference folder.
    save_path = os.path.join(
        reference_folder, os.path.basename(gff)
    )  # Gets the name of the file and writes the path.
    download_file(gff, save_path)
    with gzip.open(save_path, "rb") as f:
        file_cont = f.readlines()
    file_cont = [line.decode() for line in file_cont]
    # Gets the header.
    header = "".join(
        [line for line in file_cont if line.startswith("#") and line != "###\n"]
    )
    # Opens the rest of the lines.
    gene_loc = [line for line in file_cont if not line.startswith("#")]

    # If biotypes is not already a list, convert it to a list.
    # If biotypes is "all", convert it to a list of all biotypes in the GFF file.
    if not isinstance(biotypes, list):
        if biotypes == "all":
            # Get the last column, get to the type properties, get the value of the type properties, split all biotype values and remove duplicates.
            biotypes = list(
                set(
                    [
                        line.split("\t")[-1].split("type=")[1].split(";")[0]
                        for line in gene_loc
                    ]
                )
            )
        else:
            biotypes = [biotypes]
    # Download the id mapping file from EBI's RNAcentral, which maps each entry of the data base with
    # the corresponding entry in the database where it was obtained from.
    ref = reference_folder.split("/")[0:-1]
    ref = "/".join(ref)
    download_file(
        "http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz",
        f"{ref}/id_mapping.tsv.gz",
    )

    # Filter the RNAcentral database based on the taxonomy (9606 for human).
    if not os.path.exists(f"{reference_folder}/id_map.tsv.gz"):
        bash(
            f"zcat {ref}/id_mapping.tsv.gz | grep {tax} >{reference_folder}/id_map.tsv",
            shell=True,
        )
        os.system(f"gzip {reference_folder}/id_map.tsv")

    with gzip.open(f"{reference_folder}/id_map.tsv.gz", "rb") as gz_file:
        res = gz_file.readlines()
    res = [line.decode() for line in res if line.decode().split("\t")[3] == tax]

    # preferences = {line.split("\t")[4]:[] for line in idmapcont}
    # for pref in preferences:
    #     idmapcont_pref = list(set([line.split("\t")[1] for line in idmapcont if line.split("\t")[4] == pref]))
    #     preferences[pref].extend(idmapcont_pref)

    # For each biotype, filter the GFF file based on the prefered database.
    preferences = {"miRNA": "MIRBASE"}
    idmapcont_ndict = {}
    for pref in preferences:
        idmapcont_ndict[pref] = {
            # Keep the ID of RNAcentral and MIRBASE if it is indeed a mirbase entry.
            line.split("\t")[0]: line.split("\t")[2]
            for line in res
            if line.split("\t")[4] == pref and line.split("\t")[1] == preferences[pref]
        }
    gtf_files = {
        biotype: filter_gff(gene_loc, biotype, save_path, header, idmapcont_ndict)
        for biotype in biotypes
    }
    return gtf_files


def filter_mirbase(kegg: str, ref_file: str) -> dict[str, list[str]]:
    """
    Filter the contents of the given file based on the provided KEGG identifier.

    Args:
        kegg (str): The KEGG identifier to filter the file contents.
        ref_file (str): The reference file to be filtered.

    Returns:
        dict: A dictionary containing the filtered file contents with modifications.
    """
    # Read the contents of the reference file (the whole miRNA base database).
    with open(ref_file, "r") as r:
        fileCont = r.readlines()

    # Creates a dictionary containing only the entries with the given KEGG identifier (hsa for human) as key and the miRNA sequence as value.
    fileCont = {
        fileCont[i].rstrip(): fileCont[i + 1].rstrip()  # Removes the \n character.
        for i in range(0, len(fileCont), 2)
        if kegg in fileCont[i]
    }

    # Replaces the U with T in the miRNA sequence and saves it in a new dictionary.
    # In this new dictionary, the keys are the sequences and the values are the identifiers.
    fileCont = {
        fileCont[key].replace("U", "T"): [
            key.split(" ")[0].replace(">", ""),
            key.split(" ")[1],
        ]  # The value is both the name and the mirbase ID in a list.
        for key in fileCont
    }
    return fileCont


def get_mirna_counts(args: tuple) -> dict[str, dict[str, int | str]]:
    """
    Counts the miRNAs in the given fastq file and returns a dictionary with the sample name as key and a dictionary with the miRNAs as keys and the counts as values.

    Args:
        args (tuple): The arguments must be input in a single tuple with the following components (allows multiprocessing):
            sample_name (str): The name of the sample.
            fastq_file (str): The path to the fastq file.
            mirbaseDB (dict): A dictionary containing the miRNAs as keys and the identifiers as values.

    Returns:
        dict: A dictionary with the sample name as key and as a value, another dictionary with 2 key-value pairs. The first key is "mirna" and the value is a dictionary with the miRNA as key and the number of occurrences as value. The second key is "file" and the value is the name of the output fastq with the unassigned miRNAs.
    """

    # Unpack the arguments. The first argument is the sample name, the second is the path to the fastq file, and the third is the mirbaseDB dictionary.This allows to accept a tuple as input, making the function compatible with multiprocessing.
    sample_name, fastq_file, mirbaseDB = args

    # We use this custom function to count the miRNAs that exact match the sequence in order to reduce de computing time that requires to map and count the reads that have only partial matches to the miRNA sequences.
    # Associated reads with this functions would be not be taken into account for the mapping steps.
    # After counting the mapped reads with featureCounts, we will sum the results of both functions in order to get the total number of miRNA counts.

    lines_list = read_gzfile(
        fastq_file
    )  # Reads the fastq file and returns a list of lines.

    # After testing, the fastest way to access the data is with numpy arrays.
    # As each sequence has 4 lines of information on a fastq file,
    # the numpy array is used to split the list into chunks of 4 lines with
    # the zip function.
    gzip_cont = np.array(list(zip(*[lines_list] * 4)))

    # Creates a dictionary that unifies sequences with multiple identifiers.
    # The keys are the sequences and the values are lists of all possible identifiers.
    # mirbaseDB is the output from the filter_mirbase function.
    mirna_index = {}
    # Remember, seqs are the keys and identifiers are the values.
    for seq, mirna in mirbaseDB.items():
        # If the sequence is already in the dictionary, append the identifier to the list associated to each sequence.
        if seq in mirna_index:
            mirna_index[seq].append(mirna[0])  # The mirbase identifier is at index 0.
        # If the sequence is not in the dictionary, create a new entry with the sequence as key and the identifier as value.
        else:
            mirna_index[seq] = [mirna[0]]

    no_mirna_seqs = set()

    # Generates an empty dict ready to storage lists.
    mirna_seqs = collections.defaultdict(list)

    # This loops creates two dictionaries, mirna_seqs and no_mirna_seqs, with the sequence as key and the index of the sequence as value.
    # mirna_seqs contains the indexes of the sequences that are in mirbase, while no_mirna_seqs contains the sequences' indexes that are not in mirbase.
    for i, entry in enumerate(gzip_cont):
        # Adjust this line to match the structure of your data. Fastq files contain the sequence in the second line.
        seq = entry[1]
        # If the sequence is in mirbase, append the index of the sequence to the list associated to each sequence.
        if seq in mirna_index:
            for mirna in mirna_index[seq]:
                mirna_seqs[mirna].append(i)
        # If the sequence is not in mirbase, add its index to the list of sequences without miRNAs.
        else:
            no_mirna_seqs.add(i)

    # The following lines counts the number of occurrences of each miRNA.
    # The lenght of the list contaning the index associated to each miRNA is the number of occurrences.
    mirna_seqs_counts = {}
    for mirna in mirna_seqs:
        mirna_seqs_counts[mirna] = len(mirna_seqs[mirna])

    # We putput the sequences that are not in mirbase into a new file.
    no_mirna_seqs = [gzip_cont[i] for i in range(len(gzip_cont)) if i in no_mirna_seqs]
    outfile = f"02_trim/{sample_name}_trimmed_no_mirna.fastq.gz"

    gzip_cont = "\n".join(["\n".join(seq) for seq in no_mirna_seqs])
    with gzip.open(outfile, "wb") as f:
        f.write(gzip_cont.encode())

    # The output is a dictionary with the sample name as key and as a value, another dictionary with 2 key-value pairs.
    # The first key is "mirna" and the value is a dictionary with the miRNA as key and the number of occurrences as value.
    # The second key is "file" and the value is the name of the output fastq with the unassingned miRNAs.
    return {sample_name: {"mirna": mirna_seqs_counts, "file": outfile}}


def mirbase_sequence_assign(
    sample_dict: dict, mirbaseDB: dict, processes: int = 4
) -> tuple[dict[str, str], dict[str, int]]:
    """
    Runs get_mirna_counts in parallel using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing sample names as keys and sample data as values.
        mirbaseDB (dict): The mirBase database to be used for assigning miRNA counts and files.
        processes (int, optional): The number of processes to use for multiprocessing. Defaults to "4". If 0 is specified, use the number of samples to maximize parallelization.
    Returns:
        tuple[dict[str, str], dict[str, int]]: A tuple containing two dictionaries:
            The first dictionary maps sample names to their corresponding files.
            The second dictionary maps sample names to their corresponding miRNA counts.
    """

    if processes == 0:
        processes = len(sample_dict)

    with Pool(processes) as pool:
        sample_files = pool.map(
            get_mirna_counts,
            [
                # Creates a tuple for each sample to be used as input for the get_mirna_counts function.
                (sample_name, sample_dict[sample_name], mirbaseDB)
                for sample_name in sample_dict
            ],
        )

    # Converts sample_files and mirna_counts to standard python dictionaries.
    sample_files = dict(collections.ChainMap(*sample_files))
    mirna_counts = dict(
        collections.ChainMap(
            *[
                {sample_name: sample_files[sample_name]["mirna"]}
                for sample_name in sample_files
            ]
        )
    )
    sample_files = dict(
        collections.ChainMap(
            *[
                {sample_name: sample_files[sample_name]["file"]}
                for sample_name in sample_files
            ]
        )
    )
    return sample_files, mirna_counts


def run_aligning(args: tuple) -> dict[str, str]:
    """
    Aligns the fastq file to the index using bowtie.

    Args:
        args (tuple): The arguments must be input in a single tuple with the following components (allows multiprocessing):
            sample_name (str): The sample name.
            fastq_file (str): The fastq file path.
            index (str): The index file path.
            threads (int): The number of threads to use for alignment.
            run (bool): Whether to run the alignment.

    Returns:
        Sample dictionary with sample name as key and fastq file path as value.
    """

    # Unpack the arguments. The first argument is the sample name, the second is the fastq file path, the third is the index file path, the fourth is the number of threads to use for alignment, and the fifth the run flag. This allows to accept a tuple as input, making the function compatible with multiprocessing.
    sample_name, fastq_file, index, threads, run = args

    # After making a first direct counting using the get_mirna_counts functions in the previous step we have reduced the number of miRNAs that need to be quantify based on their mapping quality.
    # For that we use featureCounts, which requires the sample to be aligned against the reference genome.
    # Because we have previously assign a significant portion of the reads, it reduces the computing time expend on the mapping step.

    # We obtain 3 files per sample:

    # The bowtie log file, containing info about the aligment.
    logBowtie = f"00_log/{sample_name}.bowtie"
    # The sorted bam file.
    outBam = f"04_bam/{sample_name}.bam"
    # The flagstats for each bam file.
    outDedupLog = f"00_log/{sample_name}.flagstats"
    if run:
        # bowtie -f -n $mismatches_seed -e 80 -l 18 -a -m $mapping_loc --best --strata $file_genome_latest $file_reads_latest $dir/mappings.bwt\n\n";
        bash(
            f"bowtie -p {threads} -n 0 -l 18 --best --nomaqround -e 70 -k 1 -S {index} {fastq_file} 2>{logBowtie} | samtools view --threads {threads} -bS - | samtools sort --threads {threads} -o {outBam}",
            shell=True,
        )
        bash(f"samtools index {outBam}", shell=True)
        bash(f"samtools flagstat {outBam} >{outDedupLog}", shell=True)
    return {sample_name: outBam}


def align_samples(
    sample_dict: dict,
    reference: str,
    threads: int = num_threads,
    processes: int = 4,
    run: bool = False,
) -> dict[str, str]:
    """
    Runs run_aligning in parallel using multiprocessing and returns a dictionary containing the aligned samples' paths.

    Args:
        sample_dict (dict): A dictionary containing sample names as keys and sample files paths as values.
        reference (str): The reference directory.
        processes (str, optional): The number of processes for alignment. If 0 is specified, use the number of samples to maximize parallelization.
        run (bool, optional): Whether to run the alignment. Defaults to False.

    Returns:
        dict: A dictionary containing the aligned samples.
    """

    if processes == 0:
        processes = len(sample_dict)

    with Pool(processes) as pool:
        sample_dict = pool.map(
            run_aligning,
            [
                # Creates a tuple for each sample to be used as input for the run_aligning function.
                (
                    sample_name,
                    sample_dict[sample_name],
                    reference["index"],
                    threads,
                    run,
                )
                for sample_name in sample_dict
            ],
        )
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict


def get_map_quality(args: tuple) -> None:
    """
    Calculates the mapping quality for each sample.

    Args:
        args (tuple): The arguments must be input in a single tuple with the following components (allows multiprocessing):
            sample_name (str): The sample name.
            mirna_counts (dict): A dictionary containing the miRNA counts for each sample.
            run (bool): Whether to run the mapping quality calculation.
    """

    # Unpack the arguments. The first argument is the sample name, the second is the miRNA counts, and the third the run flag. This allows to accept a tuple as input, making the function compatible with multiprocessing.
    sample_name, mirna_counts, run = args

    # Gets number of reads previously assigned to miRNAs by our custom script.
    mirna_counts = sum(mirna_counts.values())

    # Reads the bowtie log file in order to read the number of mapped reads by bowtie
    # in the unsassinged miRNAs fastq file.
    with open(f"00_log/{sample_name}.bowtie", "r") as f:
        file_cont = f.readlines()

    # Gets the number of mapped reads by bowtie in the unassigned miRNAs fastq file
    # by parsing the bowtie low file and sums it with the number of counted reads.
    # It prints the ouptut in the log file.
    mapped_reads = [
        int(line.split(" ")[1]) for line in file_cont if "Reported " in line
    ][0]
    mapped_reads = mapped_reads + mirna_counts
    if run:
        log_file = f"00_log/{sample_name}.log"
        mode = "a"
        text = "\n############ MAPPING QUALITY ##############\n\n"
        write_log(log_file, text, mode)
        text = f"{sample_name} has {mapped_reads} mapped reads\n\n"
        write_log(log_file, text, mode)


def quality_mapping_samples(
    sample_dict: dict,
    mirna_counts: dict,
    processes: int = 4,
    run: bool = False,
) -> None:
    """
    Runs get_map_quality in parallel using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing sample names as keys and sample files paths as values.
        mirna_counts (dict): A dictionary containing the miRNA counts for each sample.
        processes (str, optional): The number of processes for mapping quality. If 0 is specified, use the number of samples to maximize parallelization.
        run (bool, optional): Whether to run the mapping quality calculation. Defaults to False.

    Returns:
        None
    """

    if processes == 0:
        processes = len(sample_dict)

    with Pool(processes) as pool:
        pool.map(
            get_map_quality,
            [
                # Creates a tuple for each sample to be used as input for the get_map_quality function.
                (sample_name, mirna_counts[sample_name], run)
                for sample_name in sample_dict
            ],
        )


def run_featurecount(args: tuple) -> dict[any, str]:
    """
    Runs featureCounts on the bam file and returns the output file path.

    Args:
        args (tuple): The arguments must be input in a single tuple with the following components (allows multiprocessing):
            sample_name (str): The sample name.
            bam_file (str): The bam file path.
            gff_file (str): The gff file path.
            biotype (str): The biotype to count.
            threads (int): The number of threads to use.
    """

    # Unpack the arguments. The first argument is the sample name, the second is the bam file path, the third is the gff file path, the fourth is the biotype to count, the fifth is the number of threads to use, and the sixth the run flag .This allows to accept a tuple as input, making the function compatible with multiprocessing.
    sample_name, bam_file, gff_file, biotype, threads, run = args

    out_name = f"05_counts/{sample_name}_{biotype}.counts.txt"
    if run:
        bash(
            f"featureCounts -T {threads} -t {biotype} -g Name -s 1 -O -a {gff_file} -o {out_name} {bam_file}",
            shell=True,
        )
    return {sample_name: out_name}


def quantify_biotype(
    sample_dict: dict,
    gff_file: str,
    biotype: str,
    threads: int = num_threads,
    processes: int = 4,
    run: bool = False,
) -> dict[str, str]:
    """
    Runs run_featurecount in parallel using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing sample names as keys and sample data as values.
        gff_file (str): The gff file path.
        biotype (str): The biotype to count.
        threads (int, optional): The number of threads to use. Defaults to the total number of CPU threads.
        run (bool, optional): Whether to run the evaluation. Defaults to False.
        processes (int, optional): The number of processes to use for multiprocessing. Defaults to the number of samples to maximize parallelization.

    Returns:
        A dictionary containing the quantified biotypes for each sample.
    """

    if processes == 0:
        processes = len(sample_dict)

    with Pool(processes) as pool:
        sample_dict = pool.map(
            run_featurecount,
            [
                # Creates a tuple for each sample to be used as input for the run_featurecount function.
                (
                    sample_name,
                    sample_dict[sample_name],
                    gff_file,
                    biotype,
                    threads,
                    run,
                )
                for sample_name in sample_dict
            ],
        )
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict


def quantify_mirnas(args: tuple) -> None:
    """
    This functions sums the number of assigned reads by the script and by featureCounts.

    Args:
        args (tuple): The arguments must be input in a single tuple with the following components (allows multiprocessing):
            sample_name (str): The sample name.
            mirna_counts (dict): The mirna counts dictionary.
            run (bool): Whether to run the quantification.

    Returns:
        None
    """
    # This function gets the number of assigned reads by the script and by featureCounts.

    # Unpack the arguments. The first argument is the sample name, the second is the mirna counts, and the third the run flag. This allows to accept a tuple as input, making the function compatible with multiprocessing.
    sample_name, mirna_counts, run = args

    mirna_counts = sum(mirna_counts.values())
    with open(f"05_counts/{sample_name}_miRNA.counts.txt.summary", "r") as f:
        file_cont = f.readlines()
    assigned = int(
        [line.split("\t")[1].rstrip() for line in file_cont if "Assigned" in line][0]
    )
    # Sums the number of assigned reads by the script and by featureCounts.
    assigned = assigned + mirna_counts
    if run:
        log_file = f"00_log/{sample_name}.log"
        mode = "a"
        text = "\n############ miRNA QUALITY ##############\n\n"
        write_log(log_file, text, mode)
        text = f"{sample_name} has {assigned} assigned reads\n\n"
        write_log(log_file, text, mode)


def quantify_samples(
    sample_dict: dict,
    mirna_counts: dict,
    processes: int = 4,
    run: bool = False,
) -> None:
    """
    Runs quantify_mirnas in parallel using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing sample information.
        mirna_counts (dict): A dictionary containing miRNA counts for each sample.
        processes (int, optional): The number of processes to use for multiprocessing. Defaults to 4. If 0 is specified, use the number of samples to maximize parallelization.
        run (bool, optional): Flag indicating whether the function should run. Defaults to False.

    Returns:
        None
    """

    if processes == 0:
        processes = len(sample_dict)

    with Pool(processes) as pool:
        pool.map(
            quantify_mirnas,
            [
                # Creates a tuple for each sample to be used as input for the quantify_mirnas function.
                (sample_name, mirna_counts[sample_name], run)
                for sample_name in sample_dict
            ],
        )


def concat_mirna(args: tuple) -> dict[str, str]:
    """
    This function merges the results from the featureCounts and the get_mirna_count functions.

    Args:
        args (tuple): The arguments must be input in a single tuple with the following components (allows multiprocessing):
            sample_name (str): The sample name.
            count_file (str): The count file path.
            mirna_counts (dict): A dictionary containing miRNA counts for each sample.
            use_mirbase (str): When the species we are analyzing is not in mirBase, left use_mirbase emtpy and the data from the available database will be used.
            mirbaseDB (dict): mirBase database.

    Returns:
        A dictionary containing the merged results.
    """
    # This function merges the results from the featureCounts and the get_mirna_count functions.
    # Count file is the output from featureCounts and miRNA the output from the script.
    # Remember: mirbaseDB is the output from the filter_mirbase function. The keys are the sequences and the values are the identifiers.
    # use_mirbase in this case is used as a flag. When the species we are using is not mirBase, use_mirbase is emtpy (None) and we will use the data from the available database.

    # Unpack the arguments. The first argument is the sample name, the second is the count file, the third is the mirna counts, and the fourth the use_mirbase flag. This allows to accept a tuple as input, making the function compatible with multiprocessing.
    sample_name, count_file, mirna_counts, mirbaseDB, use_mirbase = args

    print(count_file)
    with open(count_file) as f:
        read_count = f.readlines()

    # Get a list with only the count lines.
    read_count = [
        line
        for line in read_count
        if not line.startswith("#") and not line.startswith("Geneid")
    ]

    # Get only the miRNA name an the count value and create a dictionary.
    read_count = {
        line.split("\t")[0]: int(line.split("\t")[-1].rstrip()) for line in read_count
    }
    read_count = {key: read_count[key] for key in read_count if read_count[key] > 0}

    # This is the part that sums the output of each counting method in the same variable for each miRNA.
    # If use_mirbase is specified it will use mirBase.
    # If use_mirbase is empty (None), counts are identified by the available miRNA IDs for that species.
    if use_mirbase:
        mirbaseDBCounts = [
            mirbaseDB[key][0] for key in mirbaseDB
        ]  # Keep all mirbase names in a list.
        # Creates an empty dictionary where all counts will be stored.
        counts = {}
        for mirna in mirbaseDBCounts:
            counts[mirna] = 0
            # Adds the counts from featureCounts.
            if mirna in read_count:
                counts[mirna] += read_count[mirna]
            # Adds the counts from the script.
            if mirna in mirna_counts:
                counts[mirna] += mirna_counts[mirna]
        counts = {key: counts[key] for key in counts if counts[key] > 0}
    else:
        # This time storages the mirbase ID as keys and the miRNA names as values in a dictionary
        # This times uses the mirbaseID to access the output from featureCounts.
        mirbaseDBCounts = {mirbaseDB[key][1]: mirbaseDB[key][0] for key in mirbaseDB}
        counts = {}
        # mirna is the mirbaseDB ID and mirna_alt is the miRNA name.
        for mirna, mirna_alt in mirbaseDBCounts.items():
            counts[mirna_alt] = 0
            if mirna in read_count:
                # Adds the counts of the miRNAs identified with the mirbase ID.
                counts[mirna_alt] += read_count[mirna]
            if mirna_alt in mirna_counts:
                counts[mirna_alt] += mirna_counts[mirna_alt]

        # Returns only the miRNAs with counts greater than 0 combining the counts from both methods.
        counts = {key: counts[key] for key in counts if counts[key] > 0}

    # Writes the combined counts to a file.
    with open(f"05_counts/{sample_name}_miRNA_concat.txt", "w") as f:
        f.write("miRNA\tcounts\n")
        for mirna in read_count:
            f.write(f"{mirna}\t{read_count[mirna]}\n")

    # Updates the sample_dict.
    # The file path is used in the next step to assign the files to the samples.
    return {sample_name: f"05_counts/{sample_name}_miRNA_concat.txt"}


def concat_mirna_samples(
    sample_dict: dict,
    mirna_counts: dict,
    mirbaseDB: dict,
    use_mirbase: str = None,
    processes: int = 4,
) -> dict[str, str]:
    """
    Runs concat_mirna in parallel using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing sample names as keys and sample data as values.
        mirna_counts (dict): A dictionary containing miRNA counts for each sample.
        use_mirbase (str): When the species we are analyzing is not in mirBase, left use_mirbase emtpy and the data from the available database will be used.
        mirbaseDB (dict): The mirBase database to be used for assigning miRNA counts and files.
        processes (int, optional): The number of processes to use for multiprocessing. Defaults to "4". If 0 is specified, use the number of samples to maximize parallelization.

    Returns:
        A dictionary containing the concatenated miRNA counts for each sample.
    """

    if processes == 0:
        processes = len(sample_dict)

    with Pool(processes) as pool:
        sample_dict = pool.map(
            concat_mirna,
            # Creates a tuple for each sample to be used as input for the concat_mirna function.
            [
                (
                    sample_name,
                    sample_dict[sample_name],
                    mirna_counts[sample_name],
                    mirbaseDB,
                    use_mirbase,
                )
                for sample_name in sample_dict
            ],
        )
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict


def merge_count_files(
    suffix: str = None,
    folder_path: str = "05_counts/",
    sample_dict: dict = None,
    run: bool = False,
) -> None:
    """
    Merges concat count files with the given suffix and store the result in a new TSV file. If a sample_dict is provided (a dictionary containing sample names as keys and sample files paths as values), the concat count files will be merged using the sample_dict. It gives the count matrix as ouput.

    Args:
        suffix (str, optional): The file suffix to filter count files.
        folder_path (str, optional): The folder path where the count files are stored. Defaults to "05_counts/".
        sample_dict (dict, optional): A dictionary containing sample names as keys and sample files paths as values.
        run (bool, optional): Whether to run the evaluation. Defaults to False.

    Returns:
        None
    """

    if run:

        if not sample_dict:
            # Get a list of all count files in the directory
            sample_dict = {
                # Remove the suffix from the file name : folder_path + file name to create the full path.
                file.replace(suffix, ""): os.path.join(folder_path, file)
                # List the folder path.
                for file in os.listdir(folder_path)
                # Add only files with the suffix.
                if file.endswith(suffix)
            }

        # Initialize an empty DataFrame to store the merged data
        merged_data = pd.DataFrame(columns=["miRNA"])

        # Loop through each count file.
        # Sort the sample names previously.
        for sample_name in sorted(sample_dict):

            # Read the count file into a DataFrame
            df = pd.read_csv(sample_dict[sample_name], delimiter="\t")

            # Replace "-" with "_" for R compatibility.
            sample_name = sample_name.replace("-", "_")

            # Rename the columns, excluding the first column (miRNA)
            df.columns = ["miRNA"] + [sample_name for col in df.columns[1:]]

            # Merge the DataFrame with the merged_data DataFrame, using the miRNA column as the key
            merged_data = pd.merge(merged_data, df, on="miRNA", how="outer")

        # Convert all non-integer values to NaN and then replace NaN with 0
        for col in merged_data.columns[1:]:
            merged_data[col] = (
                merged_data[col]
                .apply(pd.to_numeric, errors="coerce")
                .fillna(0)
                .astype(int)
            )

        # Write the merged data to a new CSV file
        merged_data.to_csv(
            os.path.join(folder_path, "count_matrix.csv"), sep="\t", index=False
        )


def create_colData(
    groups: str,
    read_type: str,
    suffix: str = None,
    folder_path: str = "05_counts/",
    sample_dict: dict = None,
    run: bool = False,
) -> None:
    """
    Creates a colData file for DESeq2 and store it in a TSV file in the specified folder_path. To load the sample paths, a sample_dict can be provided where sample names are keys and sample files paths are values. If not, a suffix must be provided to filter the count files in the specified folder_path.

    Args:
        groups (str): A comma-separated string including the codes in the file names identifying the experiment groups.
        suffix (str, None): The file suffix to filter count files.
        read_type (str): The type of sequencing, either single-read or pair-end.
        folder_path (str, optional): The folder path where the colData file will be stored. Defaults to "05_counts/".
        sample_dict (dict, optional): A dictionary containing sample names as keys and sample files paths as values.
        run (bool, optional): Whether to run the evaluation. Defaults to False.

    Returns:
        None
    """

    if run:

        if not sample_dict:
            # Get a list of all count files in the directory
            sample_dict = {
                # Remove the suffix from the file name : folder_path + file name to create the full path.
                file.replace(suffix, ""): os.path.join(folder_path, file)
                # List the folder path.
                for file in os.listdir(folder_path)
                # Add only files with the suffix.
                if file.endswith(suffix)
            }

        # Create a list to store the data
        data = [["sample", "group", "type"]]

        # Define the regular expression pattern
        pattern = rf".*({'|'.join(groups)}).*"

        # Loop through each file.
        # Sort the sample names previously.
        for sample_name in sorted(sample_dict):
            # Match the sample name against the pattern
            match = re.match(pattern, sample_name)

            # If the sample name matches the pattern
            if match:
                sample_name = sample_name.replace("-", "_")

                # Get the group name:
                group = match.group(1)

                # Add the data to the list
                data.append([sample_name, group, read_type])

        # Create the CSV file
        with open(f"{folder_path}/colData.csv", "w") as f:
            for row in data:
                f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")
