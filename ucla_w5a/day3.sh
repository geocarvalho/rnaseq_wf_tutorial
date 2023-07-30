#!/bin/bash

cd reference/
wget ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

docker pull quay.io/biocontainers/salmon:1.10.2--hecfa306_0
docker run -v $PWD:$PWD quay.io/biocontainers/salmon:1.10.2--hecfa306_0 salmon index -t $PWD/Homo_sapiens.GRCh38.cdna.all.fa -i $PWD/human_hg19_index -k 21

cd ../
for fq in data/*_trimmed.fastq.gz; do
    base=$(basename ${fq} .fastq.gz)
    echo $base running
    docker run -v $PWD:$PWD quay.io/biocontainers/salmon:1.10.2--hecfa306_0 salmon quant \
        -i $PWD/reference/human_hg19_index \
        -l A \
        -r $PWD/$fq \
        -p 2 \
        -o $PWD/data/${base}".salmon" \
        --seqBias \
        --gcBias \
        --useVBOpt \
        --validateMappings
done

