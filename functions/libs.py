from typing import Generator, Literal


def detect_paired_single(sampleName, listFiles) -> Literal["paired", "single"]:
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


def shutil_python(output_file, input_files) -> None:
    """
    Concatenates the data from multiple input files and writes the result to the specified output file.
    Args:
        output_file (str): The path to the output file.
        input_files (List[str]): The list of paths to the input files.
    Returns:
        None
    """

    import shutil
    import subprocess

    # Concatenate the data from the input files
    with open(output_file, "wb") as out:
        for input_file in input_files:
            with open(input_file, "rb") as f_in:
                shutil.copyfileobj(f_in, out)
    # Compress the concatenated data
    subprocess.run(["gzip", output_file], check=True)


def zcat_files(output_file, input_files) -> None:
    """
    Compresses and concatenates the input files into the specified output file.

    Args:
        output_file (str): The name of the output file to be created.
        input_files (list of str): The list of input files to be concatenated and compressed.

    Returns:
        None
    """

    import os

    os.system("zcat {} >{}".format(" ".join(input_files), output_file))
    os.system("gzip {}".format(output_file))


def concatenate_files(args) -> dict[any, dict[any, any]]:
    """
    Concatenates files based on sample name and list of files and returns a dictionary
    containing sample name and corresponding files.
    This function is thought to concatenate lane-split files.

    Args:
        args (tuple): A tuple containing sample name, list of files, and a run number.

    Returns:
        dict: A dictionary with sample name as key and a nested dictionary containing
        file information as value.
    """

    import os

    sampleName, listFiles, run = args
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
        if run == "1":
            if os.path.exists(output_file + ".gz"):
                os.remove(output_file + ".gz")
            zcat_files(output_file, input_files)  # Concatenate and compress the files.
        sampleDict[strand] = output_file + ".gz"
    return {sampleName_change: sampleDict}


def mkdir(dir) -> None:
    """
    Function to create a directory
    """
    import os

    if not os.path.exists(dir):
        os.mkdir(dir)


def list_dir_files(dir, pattern="None") -> list[str]:
    """
    Function to list the files of a directory
    """

    import glob

    if pattern == "None":
        files = glob.glob(f"{dir}/*")
    else:
        files = glob.glob(f"{dir}/*{pattern}*")
    return files


def copy_files(or_file, to_file) -> None:
    """
    Function to copy files from source to dest
    """

    import shutil

    shutil.copy(or_file, to_file)


def rm_file(file) -> None:
    """
    Function to remove file if it exists
    """

    import os

    if os.path.exists(file):
        os.remove(file)


def download_file(url, filename, force=False) -> None:
    """
    Function to download files and files.gz
    """

    import os
    import requests
    import tqdm

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


def get_sample_name(file_names) -> list[any]:
    """
    Function to list the sample names of a list of fastq files
    """

    import os

    return list(
        set(
            [
                os.path.basename(file).split("_R1_")[0]
                for file in file_names
                if "_R1_" in file
            ]
        )
    )


def read_gzfile(filename) -> Generator[bytes | str, any, None]:
    # The use of a generator allows better memory management (loads the file in chunks).
    """
    Function to read a gz file
    """

    import gzip

    with gzip.open(filename, "rt") as f:
        for line in f:
            yield line.rstrip()


def write_log(logfile, text, mode) -> None:
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


def eval_fastq_file(args) -> None:
    """
    Function to evaluate a fastq file. It outputs the number of raw reads in the fastq file and the number of reads after trimming.
    Takes arguments for sample name, sample dictionary, output, adapter, threads, and run.
    """

    import numpy as np

    sample_name, sample_dict, output, adapter, threads, run = args
    import subprocess

    if run == "1":
        subprocess.run(f"fastqc {sample_dict} -o {output} -t {threads}", shell=True)
    log_file = f"00_log/{sample_name}.log"
    if adapter != "None":
        if run == "1":
            mode = "w"
            text = "############ RAW READS ##############\n\n"
            write_log(log_file, text, mode)
            mode = "a"
            filename = sample_dict
            lines = read_gzfile(filename)
            gzip_cont = np.array(list(zip(*[lines] * 4)))
            text = f"{filename} has {int(len(gzip_cont))} reads\n\n"
            write_log(log_file, text, mode)
    else:
        if run == "1":
            mode = "a"
            text = "############ TRIM READS ##############\n\n"
            write_log(log_file, text, mode)
            mode = "a"
            filename = sample_dict
            lines = read_gzfile(filename)
            gzip_cont = np.array(list(zip(*[lines] * 4)))
            text = f"{filename} has {int(len(gzip_cont))} reads\n\n"
            write_log(log_file, text, mode)


def eval_fastq_files(sample_dict, output, adapter, run, processes="sample") -> None:
    """
    Evaluate fastq files using multiprocessing for parallel processing.

    Args:
        sample_dict (dict): A dictionary containing sample information.
        output (str): The output directory for the evaluated files.
        adapter (str): The adapter to be used for evaluation.
        run (int): The run number for the evaluation.
        processes (str, optional): The number of processes to be used for multiprocessing.
            Defaults to "sample".

    Returns:
        None: This function does not return anything.
    """

    # By the default the number of processes is set to the number of samples in order to maximize the parallelization.
    if processes == "sample":
        processes = len(sample_dict)

    import multiprocessing

    num_threads = multiprocessing.cpu_count()

    with multiprocessing.Pool(processes) as pool:
        pool.map(
            eval_fastq_file,
            [
                (
                    sample_name,
                    sample_dict[sample_name],
                    output,
                    adapter,
                    num_threads,
                    run,
                )
                for sample_name in sample_dict
            ],
        )


def remove_umi_delete_adapter(fastq_file, adapter, outfile) -> int:
    """
    Remove UMIs and delete adapters from a given fastq file and save the processed data to an output file.

    Args:
    fastq_file (str): The input fastq file path.
    adapter (str): The adapter sequence to be removed. WARNING: If the adapter sequence is not correct, the function will output an empty file.
    outfile (str): The output file path to save the processed data.

    Returns:
    int: The total count of duplicated lines after processing.
    """

    import gzip
    import numpy as np

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


def run_trimming(args) -> dict[any, str]:
    """
    Trims the given fastq file using umi removal for adapter and umi removal and cutadapt for quality trimming.
    Args:
        args (tuple): A tuple containing sample_name, fastq_file, adapter, num_threads, and run.
    Returns:
        dict: A dictionary containing the sample_name as key and the trimmed fastq file as value.
    """

    import subprocess

    sample_name, fastq_file, adapter, num_threads, run = args

    # For each step, the file name is updated to indicate the process performed over the sample.
    # This file name change is registered in a variable called sample_dict. This is a dictionary with sample_name as key and fastq_file path as a value.
    # That variable is the key for the integration of all the functions in the pipeline, as it allows the correct localization of the appropriate file for each step.
    outFileUMI = f"02_trim/{sample_name}_umi.fastq.gz"
    outFileCut = f"02_trim/{sample_name}_trimmed.fastq.gz"
    if run == "1":
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
        subprocess.run(
            f"cutadapt --quiet -j {num_threads} -m 10 -M 40 -q 10 {outFileUMI} -o {outFileCut}",
            shell=True,
        )

    # Remove the intermediate files
    rm_file(outFileUMI)
    # Returns a dict containing the sample name as key and the trimmed fastq file as value.
    # This would be the basic return value for the functions in order to build the sample dict.
    # This dict relates the sample names (the identifier the functions work with) with the path of the file in order to be find by the function calls.
    # After each step, this dict is updated with the new file names for each samples as processing steps are performed over them.
    return {sample_name: outFileCut}


def trimming_files(sample_dict, adapter, run, processes="sample") -> dict[any, str]:
    """
    A function that trims files using multiprocessing and returns a dictionary containing the trimmed files.

    Args:
        sample_dict (dict): A dictionary containing the samples to be trimmed.
        adapter (str): The adapter used for trimming.
        run (int): The run number.
        processes (str, optional): The number of processes to use for multiprocessing. Defaults to "sample".

    Returns:
        dict: A dictionary containing the trimmed files.
    """

    if processes == "sample":
        processes = len(sample_dict)

    import multiprocessing
    import collections

    num_threads = multiprocessing.cpu_count()

    with multiprocessing.Pool(processes) as pool:
        sample_dict = pool.map(
            run_trimming,
            [
                (sample_name, sample_dict[sample_name], adapter, num_threads, run)
                for sample_name in sample_dict
            ],
        )

    # Collects all the generated sample_dict outputs from each call to run_trimming and creates a single dictionary where the keys are the sample name and the values the file paths for each sample. For each step, the file name is updated to indicate the performed process. In this case, teh suffix "_trimmed is added".
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict


def trimming_files_slow(sample_dict, adapter, run) -> dict[any, any]:
    """
    Trims files in the sample_dict using the specified adapter and run settings wihouth using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing sample names as keys and file paths as values.
        adapter (str): The adapter sequence for trimming.
        run (str): The run settings for trimming.

    Returns:
        dict: A dictionary containing trimmed sample names as keys and trimmed file paths as values.
    """

    from multiprocessing import cpu_count

    num_threads = cpu_count()

    # Im case the RAM is not enough, only one sample is processed at a time.
    trimmed_dict = {}
    for sample_name in sample_dict:
        trimmed_dict.update(
            run_trimming(
                args=(sample_name, sample_dict[sample_name], adapter, num_threads, run)
            )
        )

    return trimmed_dict


def convert_quality_to_numeric(quality_str) -> list[int]:
    """
    Convert ASCII quality scores to numeric values
    """

    # Convert ASCII quality scores to numeric values
    # ord obtains the ASCII value of a character.
    quality_str_num = [ord(str(char)) - 33 for char in quality_str]
    return quality_str_num


def get_fastq_stats(args):
    """
    Calculate and write quality and length statistics of FASTQ reads to log file.
    """

    sample_name, fastq, run = args
    if run == "1":
        import numpy
        import numpy as np

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
            qual = numpy.array(
                [seq[pos - 1] for seq in read_quality if len(seq) >= pos]
            )
            mean_qual = numpy.mean(qual)
            median_qual = numpy.median(qual)
            sd_qual = numpy.std(qual)
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


def get_stats_fastq_files(sample_dict, run, processes="sample") -> None:
    """
    Calculate statistics for fastq files in parallel using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing sample names as keys and
            fastq file paths as values.
        run (str): The run identifier.
        processes (str, optional): The number of processes to use. Defaults to
            "sample".

    Returns:
        None
    """

    if processes == "sample":
        processes = len(sample_dict)

    import multiprocessing

    with multiprocessing.Pool(processes) as pool:
        pool.map(
            get_fastq_stats,
            [
                (sample_name, sample_dict[sample_name], run)
                for sample_name in sample_dict
            ],
        )


def prepare_ref(fasta, ref) -> None:
    """
    Prepares reference genome index for mapping using bowtie-build.

    Args:
        fasta (str): The path to the fasta file.
        ref (str): The reference directory.

    Returns:
        None
    """

    from multiprocessing import cpu_count

    num_threads = cpu_count()

    mkdir(f"{ref}/Bowtie")
    bw_files = list_dir_files(f"{ref}/Bowtie")
    if len(bw_files) == 0:
        ########### Launch Reference ###############
        subprocess.run(f"gunzip -k {fasta}", shell=True)
        subprocess.run(
            f"bowtie-build {fasta.replace('.gz','')} {ref}/Bowtie/genome --threads {num_threads}",
            shell=True,
        )
        subprocess.run(f"rm {fasta.replace('.gz','')}", shell=True)


def filter_gff(gene_loc, biotype, save_path, header, idmapcont_ndict):
    """
    Filter gff file by biotype and save the filtered content to a new file.

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
    # Remove duplicates and format the gtf file.
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


def prepare_biotypes(reference_folder, gtf, tax, biotypes="miRNA") -> dict[str, str]:
    """
    Prepare biotypes from a given GTF file for a specific taxonomy and return the filtered GTF files.

    Args:
    - reference_folder: A string representing the directory where the GTF file will be saved.
    - gtf: A string representing the path to the input GTF file.
    - tax: A string representing the taxonomy for filtering the GTF file.
    - biotypes: A string or list of strings representing the biotypes to filter the GTF file (default is "miRNA").

    Returns:
    A dictionary containing filtered GTF files for the specified biotypes.
    """

    import gzip
    import os

    # Downloads the gtf file and saves it in the reference folder.
    save_path = os.path.join(
        reference_folder, os.path.basename(gtf)
    )  # Gets the name of the file and writes the path.
    download_file(gtf, save_path)
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
    # If biotypes is "all", convert it to a list of all biotypes in the GTF file.
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

    import subprocess

    # Filter the RNAcentral database based on the taxonomy (9606 for human).
    if not os.path.exists(f"{reference_folder}/id_map.tsv.gz"):
        subprocess.run(
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

    # For each biotype, filter the GTF file based on the prefered database.
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


def filter_mirbase(kegg, ref_file) -> dict[str, list[str]]:
    """
    Filter the contents of the given file based on the provided KEGG identifier.

    Args:
        kegg (str): The KEGG identifier to filter the file contents.
        ref_file (str): The reference file to be filtered.

    Returns:
        dict: A dictionary containing the filtered file contents with modifications.
    """
    # Read the contents of the reference file (tho whole miRNA base database).
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


def get_mirna_counts(args) -> dict[str, dict[str, any]]:
    """
    This function takes in three arguments: sample_name, fastq_file, and mirbaseDB.
    It processes the fastq_file to count the occurrences of different miRNAs based on the provided mirbaseDB.
    The function returns a dictionary containing the sample_name, mirna sequences counts, and the output file.
    """
    # We use this custom function to count the miRNAs that exact match the sequence in order to reduce de computing time that requires to map and count the reads that have only partial matches to the miRNA sequences.
    # Associated reads with this functions would be not be taken into account for the mapping steps.
    # After counting the mapped reads with featureCounts, we will sum the results of both functions in order to get the total number of miRNA counts.

    import collections
    import numpy as np

    sample_name, fastq_file, mirbaseDB = args
    lines_list = read_gzfile(
        fastq_file
    )  # Reads the fastq file and returns a list of lines.

    # After testing, the fastests way to access the data is with numpy arrays.
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
    import gzip

    gzip_cont = "\n".join(["\n".join(seq) for seq in no_mirna_seqs])
    with gzip.open(outfile, "wb") as f:
        f.write(gzip_cont.encode())

    # The output is a dictionary with the sample name as key and as a value, another dictionary with 2 key-value pairs.
    # The first key is "mirna" and the value is a dictionary with the miRNA as key and the number of occurrences as value.
    # The second key is "file" and the value is the name of the output fastq with the unassingned miRNAs.
    return {sample_name: {"mirna": mirna_seqs_counts, "file": outfile}}


def mirbase_sequence_assign(
    sample_dict, mirbaseDB, processes="sample"
) -> tuple[dict[any, any], dict[any, any]]:
    """
    Assigns miRNA counts and files to samples using multiprocessing.

    Args:
        sample_dict: A dictionary containing sample names as keys and sample data as values.
        mirbaseDB: The miRNA database to be used for assigning miRNA counts and files.
        processes: The number of processes to use for multiprocessing. Defaults to "sample".

    Returns:
        A tuple containing two dictionaries:
        - The first dictionary maps sample names to their corresponding files.
        - The second dictionary maps sample names to their corresponding miRNA counts.
    """

    if processes == "sample":
        processes = len(sample_dict)

    import multiprocessing
    import collections

    with multiprocessing.Pool(processes) as pool:
        sample_files = pool.map(
            get_mirna_counts,
            [
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


def run_aligning(args) -> dict[str, str]:
    """
    Run aligning with the given arguments using bowtie.

    Args:
        args (tuple): A tuple containing sample_name, fastq_file, index, num_threads, and run.

    Returns:
        dict: A dictionary containing the sample_name and the corresponding outBam file.
    """
    # After making a first direct counting using the get_mirna_counts functions in the previous step we have reduced the number of miRNAs that need to be quantify based on their mapping quality.
    # For that we use featureCounts, which requires the sample to be aligned against the reference genome.
    # Because we have previously assign a significant portion of the reads, it reduces the computing time expend on the mapping step.
    import subprocess

    sample_name, fastq_file, index, num_threads, run = args

    # We obtain 3 files per sample:

    # The bowtie log file, containing info about the aligment.
    logBowtie = f"00_log/{sample_name}.bowtie"
    # The sorted bam file.
    outBam = f"03_bam/{sample_name}.bam"
    # The flagstats for each bam file.
    outDedupLog = f"00_log/{sample_name}.flagstats"
    if run == "1":
        # bowtie -f -n $mismatches_seed -e 80 -l 18 -a -m $mapping_loc --best --strata $file_genome_latest $file_reads_latest $dir/mappings.bwt\n\n";
        subprocess.run(
            f"bowtie -p {num_threads} -n 0 -l 18 --best --nomaqround -e 70 -k 1 -S {index} {fastq_file} 2>{logBowtie} | samtools view --threads {num_threads} -bS - | samtools sort --threads {num_threads} -o {outBam}",
            shell=True,
        )
        subprocess.run(f"samtools index {outBam}", shell=True)
        subprocess.run(f"samtools flagstat {outBam} >{outDedupLog}", shell=True)
    return {sample_name: outBam}


def align_samples(sample_dict, reference, run, processes="sample") -> dict[str, str]:
    """
    Aligns the samples in the sample_dict to the reference using multiprocessing.

    Args:
        sample_dict (dict): A dictionary containing the samples to be aligned.
        reference (object): The reference object to which the samples will be aligned.
        run (object): The run object for the alignment process.
        processes (str, optional): The number of processes for alignment. Defaults to "sample".

    Returns:
        dict: A dictionary containing the aligned samples.
    """

    if processes == "sample":
        processes = len(sample_dict)

    import multiprocessing
    import collections

    num_threads = multiprocessing.cpu_count()

    with multiprocessing.Pool(processes) as pool:
        sample_dict = pool.map(
            run_aligning,
            [
                (
                    sample_name,
                    sample_dict[sample_name],
                    reference["index"],
                    num_threads,
                    run,
                )
                for sample_name in sample_dict
            ],
        )
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict


def get_map_quality(args) -> None:
    """
    Calculate the mapping quality based on the input arguments.

    Args:
        args (tuple): A tuple containing the sample name, miRNA counts, and run.

    Returns:
        None
    """

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
    if run == "1":
        log_file = f"00_log/{sample_name}.log"
        mode = "a"
        text = "\n############ MAPPING QUALITY ##############\n\n"
        write_log(log_file, text, mode)
        text = f"{sample_name} has {mapped_reads} mapped reads\n\n"
        write_log(log_file, text, mode)


def quality_mapping_samples(sample_dict, mirna_counts, run, processes="sample") -> None:
    """
    Map quality for each sample using multiprocessing.

    Args:
        sample_dict: A dictionary containing sample information.
        mirna_counts: A dictionary containing miRNA counts for each sample.
        run: The run information.
        processes: The number of processes to use for multiprocessing (default is "sample").

    Returns:
        None
    """

    if processes == "sample":
        processes = len(sample_dict)

    import multiprocessing

    with multiprocessing.Pool(processes) as pool:
        pool.map(
            get_map_quality,
            [
                (sample_name, mirna_counts[sample_name], run)
                for sample_name in sample_dict
            ],
        )


def run_featurecount(args) -> dict[any, str]:
    """
    Function to run featureCounts on a BAM file by calling featureCounts from subread package.

    Args:
        args (tuple): A tuple containing sample_name, bam_file, gtf_file, biotype, num_threads, and run.

    Returns:
        dict: A dictionary containing the sample_name as the key and the output file name as the value.
    """

    import subprocess

    sample_name, bam_file, gtf_file, biotype, num_threads, run = args
    out_name = f"04_counts/{sample_name}_{biotype}.counts.txt"
    if run == "1":
        subprocess.run(
            f"featureCounts -T {num_threads} -t {biotype} -g Name -s 1 -O -a {gtf_file} -o {out_name} {bam_file}",
            shell=True,
        )
    return {sample_name: out_name}


def quantify_biotype(
    sample_dict, gtf_file, biotype, run, processes="sample"
) -> dict[str, str]:
    """
    Quantifies the biotype of samples using multiprocessing.

    Args:
        sample_dict: A dictionary containing sample names as keys and sample data as values.
        gtf_file: The file path of the GTF file.
        biotype: The biotype to quantify.
        run: The run identifier.
        processes: The number of processes to use (default is "sample").

    Returns:
        A dictionary containing the quantified biotypes for each sample.
    """

    if processes == "sample":
        processes = len(sample_dict)

    import multiprocessing
    import collections

    num_trheads = multiprocessing.cpu_count()

    with multiprocessing.Pool(processes) as pool:
        sample_dict = pool.map(
            run_featurecount,
            [
                (
                    sample_name,
                    sample_dict[sample_name],
                    gtf_file,
                    biotype,
                    num_trheads,
                    run,
                )
                for sample_name in sample_dict
            ],
        )
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict


def quantify_mirnas(args) -> None:
    """
    Quantifies miRNAs based on the given arguments.

    Args:
        args (tuple): A tuple containing sample_name (str), mirna_counts (dict), and run (str).

    Returns:
        None
    """
    # This function gets the number of assigned reads by the script and by featureCounts.

    sample_name, mirna_counts, run = args
    mirna_counts = sum(mirna_counts.values())
    with open(f"04_counts/{sample_name}_miRNA.counts.txt.summary", "r") as f:
        file_cont = f.readlines()
    assigned = int(
        [line.split("\t")[1].rstrip() for line in file_cont if "Assigned" in line][0]
    )
    # Sums the number of assigned reads by the script and by featureCounts.
    assigned = assigned + mirna_counts
    if run == "1":
        log_file = f"00_log/{sample_name}.log"
        mode = "a"
        text = "\n############ miRNA QUALITY ##############\n\n"
        write_log(log_file, text, mode)
        text = f"{sample_name} has {assigned} assigned reads\n\n"
        write_log(log_file, text, mode)


def quantify_samples(sample_dict, mirna_counts, run, processes="sample") -> None:
    """
    Quantify samples using multiprocessing and a pool of processes.

    Args:
        sample_dict: A dictionary containing sample information.
        mirna_counts: A dictionary containing miRNA counts for each sample.
        run: Information about the run.
        processes (str): Number of processes to use for multiprocessing. Defaults to "sample".

    Returns:
        None
    """

    if processes == "sample":
        processes = len(sample_dict)

    import multiprocessing

    with multiprocessing.Pool(processes) as pool:
        pool.map(
            quantify_mirnas,
            [
                (sample_name, mirna_counts[sample_name], run)
                for sample_name in sample_dict
            ],
        )


def concat_mirna(args) -> dict[any, str]:
    """
    Concatenate miRNA counts from different sources and write the combined counts to a file.

    Args:
        args (tuple): A tuple containing sample_name (str), count_file (str), mirna_counts (dict), use_mirbase (str), and mirbaseDB (dict).

    Returns:
        dict: A dictionary containing the sample name as key and the file path as value.
    """
    # This function merges the results from the featureCounts and the get_mirna_count functions.
    # Count file is the output from featureCounts and miRNA the output from the script.
    # Remember: mirbaseDB si the output from the filter_mirbase function.
    # The use_mirbase is a flag to use mirBase or not for filtering the counts.

    # Unpack the arguments.
    sample_name, count_file, mirna_counts, use_mirbase, mirbaseDB = args
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
    # If use_mirbase is not 0 (true), only the counts identified with the standard miRNA name will be added (featureCounts was used with a gtf file) and Name option was specified.
    # Uf use_mirbase is 0 (false), counts from featureCounts are identified by the mirBase ID because no gtf file was used.
    if use_mirbase != "0":
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

    with open(f"04_counts/{sample_name}_miRNA_concat.txt", "w") as f:
        f.write("miRNA\tcounts\n")
        for mirna in read_count:
            f.write(f"{mirna}\t{read_count[mirna]}\n")

    return {sample_name: f"04_counts/{sample_name}_miRNA_concat.txt"}


def concat_mirna_samples(
    sample_dict, mirna_counts, use_mirbase, mirbaseDB, processes="sample"
) -> dict[any, any]:
    """
    Concatenates miRNA samples and returns a dictionary.

    Args:
        sample_dict: Dictionary containing miRNA samples
        mirna_counts: miRNA counts
        use_mirbase: Flag to use mirBase
        mirbaseDB: mirBase database
        processes: Number of processes to run in parallel (default is "sample")

    Returns:
        dict: Concatenated miRNA samples
    """

    if processes == "sample":
        processes = len(sample_dict)

    import multiprocessing
    import collections

    with multiprocessing.Pool(processes) as pool:
        sample_dict = pool.map(
            concat_mirna,
            [
                (
                    sample_name,
                    sample_dict[sample_name],
                    mirna_counts[sample_name],
                    use_mirbase,
                    mirbaseDB,
                )
                for sample_name in sample_dict
            ],
        )
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict


def merge_count_files(suffix, run, folder_path="04_counts/") -> None:
    """
    Merge count files with the given suffix and store the result in a new TSV file.
    This function is intended to be used with the "1_6_merge_count_files.py" script.
    It gives the count matrix as ouput.

    Args:
        suffix (str): The file suffix to filter count files.
        run (str): The run number.
        folder_path (str, optional): The folder path where the count files are stored. Defaults to "04_counts/".

    Returns:
        None
    """

    if run == "1":
        import os
        import pandas as pd
        import re

        # Get a list of all count files in the directory
        count_files = sorted(
            [file for file in os.listdir(folder_path) if file.endswith(suffix)]
        )

        # Initialize an empty DataFrame to store the merged data
        merged_data = pd.DataFrame(columns=["miRNA"])

        # Loop through each count file
        for file in count_files:
            # Read the count file into a DataFrame
            df = pd.read_csv(os.path.join(folder_path, file), delimiter="\t")

            # Extract the filename (excluding the file extension) to use as column header
            filename = re.sub(rf"{suffix}", "", file)
            filename = filename.replace("-", "_")

            # Rename the columns, excluding the first column (miRNA)
            df.columns = ["miRNA"] + [filename for col in df.columns[1:]]

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

        # Write the merged data to a new TSV file
        merged_data.to_csv(
            os.path.join(folder_path, "count_matrix.tsv"), sep="\t", index=False
        )


def create_colData(groups, suffix, read_type, run, folder_path="04_counts/") -> None:
    """
    Create colData file for DESeq2 based on the given groups, suffix, and read_type, and store it in a TSV file in the specified folder_path.
    Parameters:
    - groups: list of groups to match in the file names. They are used to recognise the sample group.
    - suffix: file suffix to match and remove thereafter.
    - read_type: single-read or pair-end.
    - folder_path: path of the folder where the files and TSV file will be stored (default is "04_counts/")
    Returns:
    - None
    """

    if run == "1":

        import os
        import re

        # Get all files in the folder
        files = sorted(
            [file for file in os.listdir(folder_path) if file.endswith(suffix)]
        )

        # Create a list to store the data
        data = [["sample", "group", "type"]]

        # Define the regular expression pattern
        pattern = rf".*({'|'.join(groups)}).*"

        # Loop through each file
        for file in files:
            # Match the file name against the pattern
            match = re.match(pattern, file)

            # If the file name matches the pattern
            if match:
                # Extract the filename (excluding the file extension) to use as column header
                filename = re.sub(rf"{suffix}", "", file)
                filename = filename.replace("-", "_")

                # Get the group name:
                group = match.group(1)

                # Add the data to the list
                data.append([filename, group, read_type])

        # Create the TSV file
        with open(f"{folder_path}/colData.tsv", "w") as f:
            for row in data:
                f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")
