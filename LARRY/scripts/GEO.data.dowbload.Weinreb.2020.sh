#!/bin/bash 

DATA_DOWNLOAD_DIR=/home/mvinyard/data/raw/weinreb_2020/fastqs/

GCP_ACCESSION_VALUES=$(cat /home/mvinyard/data/data_accessions/GEO.data.accession.Weinreb.2020.GoogleCloud_accessions.txt)

for ACCESSION in ${GCP_ACCESSION_VALUES}
do
  echo Now downloading FASTQ from ${ACCESSION}...
  gsutil -u broad-getzlab-mm-sclandscape -m cp -r ${ACCESSION} ${DATA_DOWNLOAD_DIR}
done
