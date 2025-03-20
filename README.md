# miRNA_Quiagen_pipe

A simple pipeline for miRNA-seq data analysis automation (from fastq to count matrix and colData for posterior DESeq2 analysis in R).
This pipeline is tailor-made to use with sequencing data generated with Qiagen's QIAseq miRNA Library Kit.

It combines a series of python scripts and calls to the linux shell for other tools such as fastqc, multiqc, bowtie and featureCounts.

For use with other library preparation kits, omit the trimming and umi removal steps and substitute them with alternative methods (cutadapt, UMItools...).

