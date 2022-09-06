#!/bin/bash

# script to perform STAR alignment to mm10

# set directories
FASTQ_DIR=/home/mvinyard/data/raw/weinreb_2020/fastqs/filtered
ALIGNED_BAM_DIR=/home/mvinyard/data/raw/weinreb_2020/bams/

# get all FASTQ file names
LIBRARIES=$(ls ${FASTQ_DIR}/*.fastq.gz)
ALL_LIBRARIES=""

for LIB in ${LIBRARIES}
do
  ALL_LIBRARIES+="${LIB}",
done

STAR --runThreadN 62 \
     --genomeDir /home/mvinyard/ref/mm10 \
     --readFilesIn ${ALL_LIBRARIES%?} \
     --readFilesCommand zcat \
     --twopassMode Basic \
     --outFileNamePrefix ${ALIGNED_BAM_DIR}\
     --outSAMtype BAM SortedByCoordinate;
