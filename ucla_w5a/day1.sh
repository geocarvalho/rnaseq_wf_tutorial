#!/bin/bash

# FASTQ files and quality control
fastqc data/*fastq.gz -o data/

# Filtering low-quality reads: The first step in sequence data analysis is usually to remove the subset of the data that has insufficient quality - keeping unreliable reads and base calls can introduce unnecessary noise in the analysis
# Manual on http://www.usadellab.org/cms/?page=trimmomatic
for fastq in $(ls data/*fastq.gz); do
    name=$(basename $fastq | cut -d. -f1)
    trimmomatic SE -phred33 ${fastq} data/${name}_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:60 -threads 4
done

# Check QC for the new FASTQ files
fastqc data/*_trimmed.fastq.gz -o data/