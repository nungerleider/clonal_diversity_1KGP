#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

ACCESSION=$1
DIR=/media/mcs/easystore1/1k_genomes
UNSORTED_BAM=${DIR}/${ACCESSION}.Aligned.out.bam
SORTED_BAM=${UNSORTED_BAM%bam}sorted.bam

echo "Downloading from the 1k Genomes project"

fasterq-dump --split-files  -e 20 -v -p $ACCESSION -O $DIR
echo "Download complete. Attempting to run STAR.."
conda activate star

STAR \
--genomeDir /home/mcs/Desktop/masked_hg38_virome \
--runThreadN 25 \
--readFilesIn ${DIR}/${ACCESSION}_1.fastq ${DIR}/${ACCESSION}_2.fastq \
--outFileNamePrefix ${DIR}/${ACCESSION}. \
--outSAMtype BAM Unsorted

echo 'STAR alignment finished.'
conda deactivate

echo "Sorting and indexing bam file.."
samtools sort $UNSORTED_BAM > $SORTED_BAM
samtools index $SORTED_BAM
samtools idxstats $SORTED_BAM > ${SORTED_BAM}.idxstats.tsv

echo 'All jobs complete'

