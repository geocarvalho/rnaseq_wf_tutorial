#!/bin/bash

# To download Ensemble human genome reference: https://useast.ensembl.org/Homo_sapiens/Info/Index
# https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/
# For gene annotations: https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/ or https://www.gencodegenes.org/human/
# Creating a STAR genome database
docker pull gvcn/star:2.7.10b
docker run -v $PWD:$PWD gvcn/star:2.7.10b STAR --help
mkdir reference; mv GRCh38_chr18.fa reference; cd reference
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
gunzip Homo_sapiens.GRCh38.109.chr.gtf.gz
gunzip gencode.v43.annotation.gtf.gz
mkdir star_database
docker run -v $PWD:$PWD gvcn/star:2.7.10b STAR --runThreadN 6 --runMode genomeGenerate --genomeFastaFiles $PWD/GRCh38_chr18.fa --sjdbGTFfile $PWD/gencode.v43.annotation.gtf --genomeDir $PWD/star_database/ --genomeSAindexNbases 12

# Aligning reads to the genome
for fastq in $(ls $PWD/data/*_trimmed.fastq.gz); do
    sample=$(basename $fastq | cut -d. -f1)
    docker run -v $PWD:$PWD gvcn/star:2.7.10b STAR --genomeDir $PWD/reference/star_database/ --readFilesIn ${fastq} --outFileNamePrefix $PWD/data/${sample}. --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 5 --outFilterMultimapNmax 1 --readFilesCommand zcat --readMapNumber 100000
done

for bam in $(ls *_trimmed.Aligned.sortedByCoord.out.bam); do docker run -v $PWD:$PWD -w $PWD biocontainers/samtools:v1.9-4-deb_cv1 samtools index $bam; done

# Tallying per-gene read counts: Use HTSeq to count the number of reads mapping to each gene 
docker pull biocontainers/htseq:v0.11.2-1-deb-py3_cv1
for bam in $(ls $PWD/data/*bam); do
    sample=$(basename $bam | cut -d. -f1)
    docker run -v $PWD:$PWD -w $PWD biocontainers/htseq:v0.11.2-1-deb-py3_cv1 htseq-count --idattr=gene_id --type=exon --mode=union --format=bam --stranded=yes ${bam} $PWD/reference/gencode.v43.annotation.gtf > $PWD/data/${sample}.pergene_counts.txt
done

# Merging all samples' counts into one big table: HTSeq generates a separate table of per-gene read counts for each sample. For statistical analyses, it is more practical to merge all these tables into one large table comprising the counts at all genes and for all samples.
# details on day2.Rmd