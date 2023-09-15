#!/bin/bash
# Align reads to the genome

# Define variables
index=$1
sampleID=$2
input_folder=$3
output_folder=$4

#STAR --genomeDir $index --genomeLoad Remove
STAR --runThreadN 4 --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn ${input_folder}/${sampleID}*R2*gz --outFileNamePrefix ${output_folder}/${sampleID} --outSAMtype BAM SortedByCoordinate --genomeLoad NoSharedMemory
STAR --genomeDir $index --genomeLoad Remove

