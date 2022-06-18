#!/bin/bash

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir /home/mvinyard/ref/mm10 \
--genomeFastaFiles /home/mvinyard/ref/mm10/refdata-gex-mm10-2020-A/fasta/genome.fa \
--sjdbGTFfile /home/mvinyard/ref/mm10/refdata-gex-mm10-2020-A/genes/genes.gtf \
--sjdbOverhang 56
