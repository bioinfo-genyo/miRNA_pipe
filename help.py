import json


def get_sample_name(file_names):
    """
    Function to list the sample names of a list of fastq files
    """
    import os

    return list(
        set([file.split("_trimmed")[0] for file in file_names if "_trimmed" in file])
    )


filenames = [
    "NR-10-POST-UDI-460_S27_trimmed.fastq.gz",
    "NR-10-PRE-UDI-451_S18_trimmed.fastq.gz",
    "NR-11-POST-UDI-454_S21_trimmed.fastq.gz",
    "NR-11-PRE-UDI-445_S12_trimmed.fastq.gz",
    "NR-17-POST-UDI-456_S23_trimmed.fastq.gz",
    "NR-17-PRE-UDI-443_S10_trimmed.fastq.gz",
    "NR-1-POST-UDI-464_S31_trimmed.fastq.gz",
    "NR-1-PRE-UDI-452_S19_trimmed.fastq.gz",
    "NR-25-POST-UDI-463_S30_trimmed.fastq.gz",
    "NR-25-PRE-UDI-442_S9_trimmed.fastq.gz",
    "NR-30-PRE-UDI-448_S15_trimmed.fastq.gz",
    "NR-3-POST-UDI-457_S24_trimmed.fastq.gz",
    "NR-3-PRE-UDI-449_S16_trimmed.fastq.gz",
    "NR-41-POST-UDI-468_S34_trimmed.fastq.gz",
    "NR-41-PRE-UDI-446_S13_trimmed.fastq.gz",
    "NR-43-PRE-UDI-436_S4_trimmed.fastq.gz",
    "NR-49-POST-UDI-439_S7_trimmed.fastq.gz",
    "NR-49-PRE-UDI-469_S35_trimmed.fastq.gz",
    "NR-4-POST-UDI-461_S28_trimmed.fastq.gz",
    "NR-4-PRE-UDI-447_S14_trimmed.fastq.gz",
    "NR-50-PRE-UDI-467_S33_trimmed.fastq.gz",
    "NR-53-POST-UDI-435_S3_trimmed.fastq.gz",
    "NR-53-PRE-UDI-471_S37_trimmed.fastq.gz",
    "NR-55-PRE-UDI-434_S2_trimmed.fastq.gz",
    "NR-59-POST-UDI-472_S38_trimmed.fastq.gz",
    "NR-59-PRE-UDI-470_S36_trimmed.fastq.gz",
    "NR-5-POST-UDI-459_S26_trimmed.fastq.gz",
    "NR-5-PRE-UDI-441_S8_trimmed.fastq.gz",
    "NR-62-PRE-UDI-465_S32_trimmed.fastq.gz",
    "NR-63-PRE-UDI-438_S6_trimmed.fastq.gz",
    "NR-65-PRE-UDI-437_S5_trimmed.fastq.gz",
    "NR-6-POST-UDI-458_S25_trimmed.fastq.gz",
    "NR-6-PRE-UDI-450_S17_trimmed.fastq.gz",
    "NR-7-POST-UDI-455_S22_trimmed.fastq.gz",
    "NR-7-PRE-UDI-453_S20_trimmed.fastq.gz",
    "NR-8-POST-UDI-462_S29_trimmed.fastq.gz",
    "NR-8-PRE-UDI-444_S11_trimmed.fastq.gz",
    "R-58-PRE-UDI-433_S1_trimmed.fastq.gz",
]


sample_names = get_sample_name(filenames)

sample_dict = {}
for sample_name in sample_names:
    fastq_file_r1 = [x for x in filenames if sample_name in x][0]
    sample_dict[sample_name] = "02_trim/" + fastq_file_r1

with open("00_log/1_2_fastq.json", "w") as jsonfile:
    json.dump(sample_dict, jsonfile, indent=4)
