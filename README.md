# miRNA_pipe

A simple pipeline for miRNA-seq data analysis automation (from fastq to count matrix and colData for posterior DESeq2 analysis in R).

It combines a series of python scripts and calls to the linux shell for other tools such as fastqc, multiqc, bowtie and featureCounts.
