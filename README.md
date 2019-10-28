# pem-insertion
pem-insertion is a pipeline built to analyze PEM-seq data(doi https://doi.org/10.1038/s41421-019-0088-8). PEM-seq is a next-generation sequencing method to quantify editing efficiency of CRISPR systems by introducing a random molecular barcode systems.
# Quick Start
## Introduction
pem-insertion consists of two main steps:
1. reads perprocessing, including deduplication according to random molecular barcode, filtering out reads with multiple constant adapter and trimming barcode and constant adapter.
2. count the number of the small insertions around the cut site
## Prerequisites
Please be sure to install:
1. python3 with pycharm(https://github.com/pysam-developers/pysam) and pandas(https://pandas.pydata.org) packages
2. BWA-mem(https://github.com/lh3/bwa)
3. FLASH(http://www.cbcb.umd.edu/software/flash)
4. seqtk(https://github.com/lh3/seqtk)
5. fastq-multx（https://github.com/brwnj/fastq-multx)
6. cutadapter（https://github.com/marcelm/cutadapt)
## Install 
git clone https://github.com/yuyang3/pem-insertion
## Running
> python pem-insertion.py -r1 R1 -r2 R2 -b BARCODE -p PRIMER -a ADAPTER -n NAME -m METAFILE -o OUTDIR -num NUMBER

METAFILE must consist of four columns seprated by "," 
for example:

translocation_name,DNA sequence of the reference,location of the cutsite in the reference,length to be specified
> perfect_1,CTGAATTAAACAGTACCATGTTCCTGGGGTGCTCAGGGAACTGCAGGGTAAAGACCCCCTAGTTCGAGAACCGCAGGGTAAAGACCCCCTCAGTTAGGGAACTGCAGGGTAAAGATATCTCCTCCCCACACTGGGGAACTGCAAGGTAGACTCCCCCCCATTTAGGGAACTGCAGGGTAAAGACCCCCCCTAGTTAGGGAACCACAGGGTAAAGATCCCCCAAGTTAGGGAATCGCAGGGTAAAGACCCCTCCCAGTAAGGGAATCACAGGGTAAAGATATCCCCTCCTCACATTAGGGAACTGCAAGGTAGATTCCCCCACAATTAGGGAACCGCCAGGTAATGACCCC,87,10
perfect_2,CTGAATTAAACAGTACCATGTTCCTGGGGTGCTCAGGGAACTGCAGGGTAAAGACCCCCTAGTTCGAGAACCGCAGGGTAAAGACCCCCTCAGTTAGGGAACTGCAGGGTAAAGATATCTCCTCCCCACACTGGGGAACTGCAAGGTAGACTCCCCCCCATTTAGGGAACTGCAGGGTAAAGACCCCCCCTAGTTAGGGAACCACAGGGTAAAGATCCCCCAAGTTAGGGAATCGCAGGGTAAAGACCCCTCCCAGTAAGGGAATCACAGGGTAAAGATATCCCCTCCTCACATTAGGGAACTGCAAGGTAGATTCCCCCACAATTAGGGAACCGCCAGGTAATGACCCC,88,10

# Acknowledgments
This pipeline is developed based on superQ(https://github.com/liumz93/superQ)
