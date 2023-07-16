#!/bin/bash
# Video https://www.youtube.com/watch?v=lG11JjovJHE
# Download data from https://www.youtube.com/redirect?event=video_description&redir_token=QUFFLUhqbTZoZ3NDOE5qU3dPSUhITFRfSWt6djZxeUNEUXxBQ3Jtc0ttY2FjSVVDbGJoY1ZURUVXNWgtbU1JWDF3YmlYNms5LXA2UUpjVE84d1BfR3dFNVN0UnpfSGZjQ0JkRUJ6TGNUSklfRGJrcnNuNFFMcEo1YjgzbDFFZ0pwM1pUVTMtR080LUJFc2tIejZYd3JNRElkYw&q=https%3A%2F%2Fdrive.google.com%2Ffile%2Fd%2F1DGHjbhcRy_zTm6H9C_AUpkzBML-JhtA3%2Fview%3Fusp%3Dsharing&v=lG11JjovJHE

# STEP 1: Run fastqc
# conda install -c bioconda fastqc
fastqc data/demo.fastq.gz -o data/

# run trimmomatic to trim reads with poor quality
# conda install -c bioconda trimmomatic
trimmomatic SE -threads 4 data/demo.fastq.gz data/demo_trimmed.fastq.gz TRAILING:10 -phred33
echo "Trimmomatic finished running!"
fastqc data/demo_trimmed.fastq.gz -o data/

# STEP 2: Run HISAT2
mkdir HISAT2
# get the genome indices http://daehwankimlab.github.io/hisat2/download/
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzf grch38_genome.tar.gz
mv grch38/ HISAT2/

# run alignment
# conda install -c bioconda hisat2
# docker pull quay.io/biocontainers/hisat2:2.2.1--h87f3376_5
# docker pull biocontainers/samtools:v1.9-4-deb_cv1
docker run -v $PWD:$PWD quay.io/biocontainers/hisat2:2.2.1--h87f3376_5 hisat2 \
    -q --rna-strandness R \
    -x $PWD/HISAT2/grch38/genome \
    -U $PWD/data/demo_trimmed.fastq.gz > $PWD/HISAT2/demo_trimmed.sam

docker run -v $PWD:$PWD docker.io/biocontainers/samtools:v1.9-4-deb_cv1 samtools sort $PWD/HISAT2/demo_trimmed.sam -o $PWD/HISAT2/demo_trimmed.bam
rm $PWD/HISAT2/demo_trimmed.sam 
echo "HISAT2 finished running!"


# STEP 3: Run featureCounts - Quantification
# get gtf
# docker pull quay.io/biocontainers/subread:2.0.6--he4a0461_0
# Where to find it: https://useast.ensembl.org/Homo_sapiens/Info/Index?db=core
wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
gunzip Homo_sapiens.GRCh38.106.gtf.gz
mv Homo_sapiens.GRCh38.106.gtf HISAT2/grch38/Homo_sapiens.GRCh38.106.gtf
docker run -v $PWD:$PWD quay.io/biocontainers/subread:2.0.6--he4a0461_0 featureCounts \
    -S 2 -a $PWD/HISAT2/grch38/Homo_sapiens.GRCh38.106.gtf \
    -o $PWD/quants/demo_featurecounts.txt \
    $PWD/HISAT2/demo_trimmed.bam
echo "featureCounts finished running!"