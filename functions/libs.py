def mkdir(dir):
    """
        Function to create a directory 
    """
    import os
    if not os.path.exists(dir):
        os.mkdir(dir)

def list_dir_files(dir,pattern = "None"):
    """
        Function to list the files of a directory
    """
    import glob
    if pattern == "None":
        files = glob.glob(f"{dir}/*")
    else:
        files = glob.glob(f"{dir}/*{pattern}*")
    return files

def copy_files(or_file,to_file):
    """
        Function to copy files from source to dest
    """
    import shutil
    shutil.copy(or_file,to_file)

def rm_file(file):
    """
        Function to remove file if it exists
    """
    import os
    if os.path.exists(file):
        os.remove(file)

def download_file(url, filename, force = False):
    """
        Function to download files and files.gz
    """
    import os
    import requests
    import tqdm
    if (not os.path.exists(filename) and not os.path.exists(filename.replace(".gz","")) and not os.path.exists(filename+".gz")) or force:
        rm_file(filename)
        rm_file(filename.replace(".gz",""))
        with open(filename, 'wb') as f:
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                total = int(r.headers.get('content-length', 0))

            # tqdm has many interesting parameters. Feel free to experiment!
                tqdm_params = {
                    'desc': url,
                    'total': total,
                    'miniters': 1,
                    'unit': 'B',
                    'unit_scale': True,
                    'unit_divisor': 1024,
                }
                with tqdm.tqdm(**tqdm_params) as pb:
                    for chunk in r.iter_content(chunk_size=8192):
                        pb.update(len(chunk))
                        f.write(chunk)

def get_sample_name(file_names):
    """
        Function to list the sample names of a list of fastq files
    """
    import os
    return(list(set([os.path.basename(file).split("_R1_")[0] for file in file_names if "_R1_" in file])))

def read_gzfile(filename):
    """
        Function to read a gz file
    """
    import gzip
    with gzip.open(filename, "rt") as f:
        for line in f:
            yield line.rstrip()

def write_log(logfile,text,mode):
    with open(logfile,mode) as write_file:
        write_file.write(text)


def eval_fastq_file(args):
    import numpy as np
    sample_name,sample_dict,output,adapter,threads,run = args
    import subprocess
    if run == "1":
        subprocess.run(f"fastqc {sample_dict} -o {output} -t {threads}",shell=True)
    log_file = f"00_log/{sample_name}.log"
    if adapter != "None":
        if run == "1":
            mode = "w"
            text = "############ RAW READS ##############\n\n"
            write_log(log_file,text,mode)
            mode = "a"
            filename = sample_dict
            lines = read_gzfile(filename)
            gzip_cont = np.array(list(zip(*[lines] * 4)))
            text = f"{filename} has {int(len(gzip_cont))} reads\n\n"
            write_log(log_file,text,mode)
    else:
        if run == "1":
            mode = "a"
            text = "############ TRIM READS ##############\n\n"
            write_log(log_file,text,mode)
            mode = "a"
            filename = sample_dict
            lines = read_gzfile(filename)
            gzip_cont = np.array(list(zip(*[lines] * 4)))
            text = f"{filename} has {int(len(gzip_cont))} reads\n\n"
            write_log(log_file,text,mode)

def eval_fastq_files(sample_dict,output,adapter,run):
    import multiprocessing
    with multiprocessing.Pool(len(sample_dict)) as pool:
        pool.map(eval_fastq_file,[(sample_name,sample_dict[sample_name],output,adapter,8,run) for sample_name in sample_dict])


def remove_umi_delete_adapter(fastq_file,adapter,outfile):
    import gzip
    import numpy as np
    lines_list = read_gzfile(fastq_file)
    gzip_cont = np.array(list(zip(*[lines_list] * 4)))
    unique_elements = set()
    filtered_lines = []
    duplicated_lines = []
    for sublist in gzip_cont:
        seq, quality = sublist[1], sublist[3]
        if seq not in unique_elements:
            if adapter in seq:
                unique_elements.add(seq)
                idx_pos = seq.rfind(adapter)
                seq_wo_adapter = seq[:idx_pos]
                quality_wo_adapter = quality[:idx_pos]
                sublist[1] = seq_wo_adapter
                sublist[3] = quality_wo_adapter
                filtered_lines.append(sublist)
        else:
            duplicated_lines.append(seq)
    unique_dup = list(set(duplicated_lines))
    duplicated_lines = len(duplicated_lines) + len(unique_dup)
    gzip_cont = "\n".join(["\n".join(seq) for seq in filtered_lines])
    with gzip.open(outfile, "wb") as f:
        f.write(gzip_cont.encode())
    return(duplicated_lines)

def run_trimming(args):
    import subprocess
    sample_name, fastq_file, adapter, num_threads, run = args
    outFileUMI = f"02_trim/{sample_name}_umi.fastq.gz"
    outFileCut = f"02_trim/{sample_name}_trimmed.fastq.gz"
    if run == "1":
        print("Running umi removal {}".format(fastq_file))
        duplicated_lines = remove_umi_delete_adapter(fastq_file,adapter,outFileUMI)
        mode = "a"
        text = "############ DUPS READS ##############\n\n"
        log_file = f"00_log/{sample_name}.log"
        write_log(log_file,text,mode)
        mode = "a"
        text = f"{sample_name} has {duplicated_lines} duplicated reads\n\n"
        write_log(log_file,text,mode)
        print("Running cutadapt {}".format(fastq_file))
        subprocess.run(f"cutadapt --quiet -j {num_threads} -m 10 -M 40 -q 10 {outFileUMI} -o {outFileCut}",shell=True)
        
    rm_file(outFileUMI)
    return({sample_name:outFileCut})

def trimming_files(sample_dict,adapter,run):
    import multiprocessing
    import collections
    with multiprocessing.Pool(len(sample_dict)) as pool:
        sample_dict = pool.map(run_trimming,[(sample_name,sample_dict[sample_name],adapter,8,run) for sample_name in sample_dict])
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return(sample_dict)

def convert_quality_to_numeric(quality_str):
    # Convert ASCII quality scores to numeric values
    quality_str_num = [ord(str(char)) - 33 for char in quality_str]
    return(quality_str_num)

def get_fastq_stats(args):
    sample_name,fastq,run = args
    if run == "1":
        import numpy
        import numpy as np
        lines = read_gzfile(fastq)
        gzip_cont = np.array(list(zip(*[lines] * 4)))
        quali_seq = [line[3].rstrip() for line in gzip_cont]
        read_quality = [convert_quality_to_numeric(seq) for seq in quali_seq]
        read_lengths = [len(qual) for qual in read_quality]
        max_i = max(read_lengths)
        quality_stats = {}
        lengths_stats = {}
        for pos in range(1, max_i + 1):
            qual = numpy.array([seq[pos-1] for seq in read_quality if len(seq) >= pos])
            mean_qual = numpy.mean(qual)
            median_qual = numpy.median(qual)
            sd_qual = numpy.std(qual)
            quality_stats[str(pos)] = {
                "nreads": len(qual),
                "median": median_qual,
                "mean": mean_qual,
                "sd": sd_qual
                }
            lengths_stats[str(pos)] = read_lengths.count(pos)
    
        logfile = f"00_log/{sample_name}.log"
        mode = "a"
        ### quality stats
        text = "############ QUALITY READS ##############\n\n"
        write_log(logfile,text,mode)
        text = "bp\tmean\n"
        write_log(logfile,text,mode)
        for i in quality_stats:
            text = f"{i}\t{quality_stats[i]['mean']}\n"
            write_log(logfile,text,mode)
    
        ### lengths stats
        text = "\n############ LEN READS ##############\n\n"
        write_log(logfile,text,mode)
        text = "len\tdensity\n"
        write_log(logfile,text,mode)
        for i in lengths_stats:
            text = f"{i}\t{lengths_stats[i]}\n"
            write_log(logfile,text,mode)


def get_stats_fastq_files(sample_dict,run):
    import multiprocessing
    with multiprocessing.Pool(len(sample_dict)) as pool:
        pool.map(get_fastq_stats,[(sample_name,sample_dict[sample_name],run) for sample_name in sample_dict])


def prepare_ref(fasta,ref):
    import subprocess
    mkdir(f"{ref}/Bowtie")
    bw_files = list_dir_files(f"{ref}/Bowtie")
    if len(bw_files) == 0:
        ########### Launch Reference ###############
        subprocess.run(f"gunzip {fasta}",shell=True)
        subprocess.run(f"bowtie-build {fasta.replace('.gz','')} {ref}/Bowtie/genome --threads 10",shell=True)
        subprocess.run(f"gzip {fasta.replace('.gz','')}",shell=True)


def filter_gff(gene_loc, biotype,save_path,header,idmapcont_ndict):
    # Filter biotype
    gene_loc_biotype = [line for line in gene_loc if line.split("\t")[-1].split("type=")[1].split(";")[0] == biotype]
    # List to dictionary
    gene_loc_biotype_dict = {line.split("\t")[-1].split("Name=")[1].split("_")[0]:line for line in gene_loc_biotype}
    # Filter ID
    idmapcont_ndict = idmapcont_ndict[biotype]
    idmapcont_ndict = {key:idmapcont_ndict[key] for key in idmapcont_ndict if key in gene_loc_biotype_dict}
    gene_loc_biotype_dict = {key:gene_loc_biotype_dict[key] for key in gene_loc_biotype_dict if key in idmapcont_ndict}
    # Remove duplicates
    gene_loc_biotype = []
    for key in gene_loc_biotype_dict:
        line = gene_loc_biotype_dict[key].split("\t")
        line[2] = biotype
        line[-1] = "Name={}".format(key)
        line = "\t".join(line)
        if line not in gene_loc_biotype:
            gene_loc_biotype.append(line)
    
    dict_map = []
    for line in gene_loc_biotype:
        line = line.split("\t")
        name = line[-1].replace("Name=","")
        nname = idmapcont_ndict[name]
        line[-1] = "Name={}".format(nname)
        line = "\t".join(line)
        dict_map.append(line)

    gtf_cont = "\n".join(dict_map)
    gtf_cont = f"{header}{gtf_cont}"
    outfile = save_path.replace(".gff3.gz","_{}.gff3".format(biotype))
    with open(outfile,"w") as out:
        out.write(gtf_cont)
    return(outfile)

def prepare_biotypes(reference_folder,gtf,tax,biotypes = "miRNA"):
    import gzip
    import os
    save_path = os.path.join(reference_folder,os.path.basename(gtf))
    download_file(gtf,save_path)
    with gzip.open(save_path,"rb") as f:
        file_cont = f.readlines()
    file_cont = [line.decode() for line in file_cont]
    header = "".join([line for line in file_cont if line.startswith("#") and line != "###\n"])
    gene_loc = [line for line in file_cont if not line.startswith("#")]
    if not isinstance(biotypes,list):
        if biotypes == "all":
            biotypes = list(set([line.split("\t")[-1].split("type=")[1].split(";")[0] for line in gene_loc]))
        else:
            biotypes = [biotypes]
    


    ref = reference_folder.split("/")[0:-1]
    ref = "/".join(ref)
    download_file("http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/id_mapping.tsv.gz",f"{ref}/id_mapping.tsv.gz")
    
    import subprocess
    if not os.path.exists(f"{reference_folder}/id_map.tsv.gz"):
        subprocess.run(f"zcat {ref}/id_mapping.tsv.gz | grep {tax} >{reference_folder}/id_map.tsv",shell=True)
        os.system(f"gzip {reference_folder}/id_map.tsv")
    
    with gzip.open(f"{reference_folder}/id_map.tsv.gz","rb") as gz_file:
        res = gz_file.readlines()
    res = [line.decode() for line in res if line.decode().split("\t")[3] == tax]

    # preferences = {line.split("\t")[4]:[] for line in idmapcont}
    # for pref in preferences:
    #     idmapcont_pref = list(set([line.split("\t")[1] for line in idmapcont if line.split("\t")[4] == pref]))
    #     preferences[pref].extend(idmapcont_pref)
    
    preferences = {"miRNA":"MIRBASE"}
    idmapcont_ndict = {}
    for pref in preferences:
        idmapcont_ndict[pref] = {line.split("\t")[0]:line.split("\t")[2] for line in res if line.split("\t")[4] == pref and line.split("\t")[1] == preferences[pref]}
    gtf_files = {biotype:filter_gff(gene_loc,biotype,save_path,header,idmapcont_ndict) for biotype in biotypes}
    return(gtf_files)

def filter_mirbase(kegg,ref_file):
    with open(ref_file,"r") as r:
        fileCont = r.readlines()
    fileCont = {fileCont[i].rstrip():fileCont[i+1].rstrip() for i in range(0,len(fileCont),2) if kegg in fileCont[i]}
    fileCont = {fileCont[key].replace("U","T"):[key.split(" ")[0].replace(">",""),key.split(" ")[1]] for key in fileCont}
    return(fileCont)

def get_mirna_counts(args):
    import collections
    import numpy as np
    sample_name, fastq_file, mirbaseDB = args
    lines_list = read_gzfile(fastq_file)
    gzip_cont = np.array(list(zip(*[lines_list] * 4)))

    mirna_index = {}
    for seq, mirna in mirbaseDB.items():
        if seq in mirna_index:
            mirna_index[seq].append(mirna[0])
        else:
            mirna_index[seq] = [mirna[0]]
    
    no_mirna_seqs = set()
    mirna_seqs = collections.defaultdict(list)

    for i, entry in enumerate(gzip_cont):
        seq = entry[1]  # Adjust this line to match the structure of your data
        if seq in mirna_index:
            for mirna in mirna_index[seq]:
                mirna_seqs[mirna].append(i)
        else:
            no_mirna_seqs.add(i)
    
    mirna_seqs_counts = {}
    for mirna in mirna_seqs:
        mirna_seqs_counts[mirna] = len(mirna_seqs[mirna])

    no_mirna_seqs = [gzip_cont[i] for i in range(len(gzip_cont)) if i in no_mirna_seqs]
    outfile = f"02_trim/{sample_name}_trimmed_no_mirna.fastq.gz"
    import gzip
    gzip_cont = "\n".join(["\n".join(seq) for seq in no_mirna_seqs])
    with gzip.open(outfile, "wb") as f:
        f.write(gzip_cont.encode())

    return({sample_name:{"mirna":mirna_seqs_counts,"file":outfile}})

def mirbase_sequence_assign(sample_dict, mirbaseDB):
    import multiprocessing
    import collections
    with multiprocessing.Pool(len(sample_dict)) as pool:
        sample_files = pool.map(get_mirna_counts,[(sample_name,sample_dict[sample_name],mirbaseDB) for sample_name in sample_dict])
    
    sample_files = dict(collections.ChainMap(*sample_files))
    mirna_counts = dict(collections.ChainMap(*[{sample_name:sample_files[sample_name]["mirna"]} for sample_name in sample_files]))
    sample_files = dict(collections.ChainMap(*[{sample_name:sample_files[sample_name]["file"]} for sample_name in sample_files]))
    return sample_files, mirna_counts


def run_aligning(args):
    import subprocess
    sample_name, fastq_file,index,num_threads, run = args
    logBowtie = f"00_log/{sample_name}.bowtie"
    outBam = f"03_bam/{sample_name}.bam"
    outDedupLog = f"00_log/{sample_name}.flagstats"
    if run == "1":
        # bowtie -f -n $mismatches_seed -e 80 -l 18 -a -m $mapping_loc --best --strata $file_genome_latest $file_reads_latest $dir/mappings.bwt\n\n";
        subprocess.run(f"bowtie -p {num_threads} -n 0 -l 18 --best --nomaqround -e 70 -k 1 -S {index} {fastq_file} 2>{logBowtie} | samtools view --threads {num_threads} -bS - | samtools sort --threads {num_threads} -o {outBam}",shell=True)
        subprocess.run(f"samtools index {outBam}",shell=True)
        subprocess.run(f"samtools flagstat {outBam} >{outDedupLog}",shell=True)
    return({sample_name:outBam})

def align_samples(sample_dict, reference, run):
    import multiprocessing
    import collections
    with multiprocessing.Pool(len(sample_dict)) as pool:
        sample_dict = pool.map(run_aligning,[(sample_name,sample_dict[sample_name],reference["index"],8,run) for sample_name in sample_dict])
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return(sample_dict)

def get_map_quality(args):
    sample_name, mirna_counts, run = args
    mirna_counts = sum(mirna_counts.values())
    
    with open(f"00_log/{sample_name}.bowtie","r") as f:
        file_cont = f.readlines()
    
    mapped_reads = [int(line.split(" ")[1]) for line in file_cont if "Reported " in line][0]
    mapped_reads = mapped_reads + mirna_counts
    if run == "1":
        log_file = f"00_log/{sample_name}.log"
        mode = "a"
        text = "\n############ MAPPING QUALITY ##############\n\n"
        write_log(log_file, text, mode)
        text = f"{sample_name} has {mapped_reads} mapped reads\n\n"
        write_log(log_file, text, mode)

def quality_mapping_samples(sample_dict, mirna_counts, run):
    import multiprocessing
    with multiprocessing.Pool(len(sample_dict)) as pool:
        pool.map(get_map_quality,[(sample_name,mirna_counts[sample_name], run) for sample_name in sample_dict])

def run_featurecount(args):
    import subprocess
    sample_name, bam_file, gtf_file,biotype, num_threads, run = args
    out_name = f"04_counts/{sample_name}_{biotype}.counts.txt"
    if run == "1":
        subprocess.run(f"featureCounts -T {num_threads} -t {biotype} -g Name -s 1 -O -a {gtf_file} -o {out_name} {bam_file}",shell=True)
    return({sample_name:out_name})


def quantify_biotype(sample_dict, gtf_file, biotype, run):
    import multiprocessing
    import collections
    with multiprocessing.Pool(len(sample_dict)) as pool:
        sample_dict = pool.map(run_featurecount,[(sample_name,sample_dict[sample_name],gtf_file, biotype, 8, run) for sample_name in sample_dict])
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict


def quantify_mirnas(args):
    sample_name, mirna_counts, run = args
    mirna_counts = sum(mirna_counts.values())
    with open(f"04_counts/{sample_name}_miRNA.counts.txt.summary","r") as f:
        file_cont = f.readlines()
    assigned = int([line.split("\t")[1].rstrip() for line in file_cont if "Assigned" in line][0])
    assigned = assigned + mirna_counts
    if run == "1":
        log_file = f"00_log/{sample_name}.log"
        mode = "a"
        text = "\n############ miRNA QUALITY ##############\n\n"
        write_log(log_file, text, mode)
        text = f"{sample_name} has {assigned} assigned reads\n\n"
        write_log(log_file, text, mode)

def quantify_samples(sample_dict, mirna_counts, run):
    import multiprocessing
    with multiprocessing.Pool(len(sample_dict)) as pool:
        pool.map(quantify_mirnas,[(sample_name,mirna_counts[sample_name], run) for sample_name in sample_dict])




def concat_mirna(args):
    sample_name, count_file, mirna_counts, use_mirbase, mirbaseDB = args
    print(count_file)
    with open(count_file) as f:
        read_count = f.readlines()
    
    read_count = [line for line in read_count if not line.startswith("#") and not line.startswith("Geneid")]
    read_count = {line.split("\t")[0]:int(line.split("\t")[-1].rstrip()) for line in read_count}
    read_count = {key:read_count[key] for key in read_count if read_count[key] > 0}

    if use_mirbase != "0":
        mirbaseDBCounts = [mirbaseDB[key][0] for key in mirbaseDB]
        counts = {}
        for mirna in mirbaseDBCounts:
            counts[mirna] = 0
            if mirna in read_count:
                counts[mirna] += read_count[mirna]
            if mirna in mirna_counts:
                counts[mirna] += mirna_counts[mirna]
        counts = {key:counts[key] for key in counts if counts[key] > 0}
    else:
        mirbaseDBCounts = {mirbaseDB[key][1]:mirbaseDB[key][0] for key in mirbaseDB}
        counts = {}
        for mirna, mirna_alt in mirbaseDBCounts.items():
            counts[mirna_alt] = 0
            if mirna in read_count:
                counts[mirna_alt] += read_count[mirna]
            if mirna_alt in mirna_counts:
                counts[mirna_alt] += mirna_counts[mirna_alt]
        counts = {key:counts[key] for key in counts if counts[key] > 0}

    with open(f"04_counts/{sample_name}_miRNA_concat.txt","w") as f:
        f.write("miRNA\tcounts\n")
        for mirna in read_count:
            f.write(f"{mirna}\t{read_count[mirna]}\n")

    return({sample_name:f"04_counts/{sample_name}_miRNA_concat.txt"})

def concat_mirna_samples(sample_dict, mirna_counts, use_mirbase, mirbaseDB):
    import multiprocessing
    import collections
    with multiprocessing.Pool(len(sample_dict)) as pool:
        sample_dict = pool.map(concat_mirna,[(sample_name,sample_dict[sample_name],mirna_counts[sample_name], use_mirbase, mirbaseDB) for sample_name in sample_dict])
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return sample_dict