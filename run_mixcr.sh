#!/bin/bash

###
### Take paired end fastq.gz files as input and run the entire MiXCR pipeline
### resulting in a final output file of SAMPLE_ID.clones.txt.
###
### Convert SAMPLE_ID.clones.txt to vdjtools compatible format 
### via the "vdjtools Convert" accessory function
###

mixcr align -p rna-seq -s hsa -OallowPartialAlignments=true ${1} ${1%1.fastq.gz}2.fastq.gz ${1%fastq.gz}.alignments.vdjca
mixcr assemblePartial ${1%fastq.gz}.alignments.vdjca ${1%fastq.gz}.alignments_rescued_1.vdjca
mixcr assemblePartial ${1%fastq.gz}.alignments_rescued_1.vdjca ${1%fastq.gz}.alignments_rescued_2.vdjca
mixcr extend ${1%fastq.gz}.alignments_rescued_2.vdjca ${1%fastq.gz}.alignments_rescued_2_extended.vdjca
mixcr assemble ${1%fastq.gz}.alignments_rescued_2_extended.vdjca ${1%fastq.gz}.clones.clns
mixcr exportClones ${1%fastq.gz}.clones.clns ${1%fastq.gz}.clones.txt

java -Xmx16G -jar vdjtools.jar Convert -S MiXCR ${1%fastq.gz}.clones.txt /media/mcs/easystore2/1000_genomes_RNA/'''
