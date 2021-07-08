# LARRY data preprocessing and analysis

## Data download

Processed data downloaded from the Klein Lab's [**paper-data repository** (commit d4528a4)](https://github.com/AllonKleinLab/paper-data/blob/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation/README.md)

Using **`gsutils`**, I first downloaded the FASTQ files uploaded to the SRA. These files had already been processed in part by the authors using [**this script**](https://github.com/AllonKleinLab/LARRY/blob/master/LARRY_sorting_and_filtering.py) where the input files are those of the inDrops pipeline. 

To download all of the google-stored data, I first wrote a web-scraping function (shown in this [**notebook**](https://github.com/mvinyard/data/blob/LARRY/LARRY/notebooks/GSM4185642_get_GCP_accessions.ipynb)) to search for each SRR value and find the corresponding GCP storage path. , concatentate them into a list, and save to a [file](https://) to be picked up by the following script:

>Google Cloud Accessions saved to disk at: /home/mvinyard/data/data_accessions/GEO.data.accession.Weinreb.2020.GoogleCloud_accessions.txt

The script detailing the **`gsutils`** download:

```BASH
#!/bin/bash 
# GEO.data.dowbload.Weinreb.2020.sh

DATA_DOWNLOAD_DIR=/home/mvinyard/data/raw/weinreb_2020/fastqs/

GCP_ACCESSION_VALUES=$(cat /home/mvinyard/data/data_accessions/GEO.data.accession.Weinreb.2020.GoogleCloud_accessions.txt)

for ACCESSION in ${GCP_ACCESSION_VALUES}
do
  echo Now downloading FASTQ from ${ACCESSION}...
  gsutil -u REDACTED_PROJECT_ID -m cp -r ${ACCESSION} ${DATA_DOWNLOAD_DIR}
done
```

Run using:
```BASH
bash GEO.data.dowbload.Weinreb.2020.sh
```

## FASTQ filtering

To streamline preprocessing and ensure an accurate reproduction of the author's results, I decided to filter the FASTQ wherein only reads with cell barcodes present in the final dataset are retained in the filtered raw reads. 

This [**notebook**](https://github.com/mvinyard/data/blob/LARRY/LARRY/notebooks/Filter_LARRY_FASTQs.ipynb) shows an example of a single library (equates to a single FASTQ file) being filtered. All libraries were filtered with the script, [filter.fastqs.Weinreb.2020.py](https://github.com/mvinyard/data/blob/LARRY/LARRY/scripts/filter.fastqs.Weinreb.2020.py) and run as follows:

```BASH
python filter.fastqs.Weinreb.2020.py
```

I created a [**log file**](https://github.com/mvinyard/data/blob/LARRY/LARRY/files/logs/filter.fastqs.Weinreb.2020.debug.log) to document this filtering process. 

## Read alignment

I first downloaded the mm10 reference files from [10x Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) using:
```BASH 
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"
```

I generated a STAR-compatible reference genome and annotation running STAR in [**genomeGenerate**](https://github.com/mvinyard/data/blob/LARRY/LARRY/scripts/STAR.genome_generate.mm10.sh) mode:
```BASH
#!/bin/bash

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir /home/mvinyard/ref/mm10 \
--genomeFastaFiles /home/mvinyard/ref/mm10/refdata-gex-mm10-2020-A/fasta/genome.fa \
--sjdbGTFfile /home/mvinyard/ref/mm10/refdata-gex-mm10-2020-A/genes/genes.gtf \
--sjdbOverhang 56
```

I then aligned all libraries using **[STAR](https://github.com/mvinyard/data/blob/LARRY/LARRY/scripts/STAR.align.mm10.LARRY.sh)**:

***Quick note***: run this command before beginning. Oftentimes, when using many cores (e.g., 32), I've experienced errors as a result of this parameter being set too low. In essence, enlarging this parameter enables your machine to manage more open files at a single time. 
```BASH
ulimit -n 147456 # (default: 1024 * 144)
```
The STAR script for alignment is as follows:
```BASH 
#!/bin/bash
# script to perform STAR alignment to mm10

# set directories
FASTQ_DIR=/home/mvinyard/data/raw/weinreb_2020/fastqs/filtered
ALIGNED_BAM_DIR=/home/mvinyard/data/raw/weinreb_2020/bams/

# get all FASTQ file names
LIBRARIES=$(ls ${FASTQ_DIR}/*.fastq)
ALL_LIBRARIES=""

for LIB in ${LIBRARIES}
do
  ALL_LIBRARIES+="${LIB}",
done

STAR --runThreadN 30 \
     --genomeDir /home/mvinyard/ref/mm10 \
     --readFilesIn ${ALL_LIBRARIES%?} \
     --twopassMode Basic \
     --outFileNamePrefix ${ALIGNED_BAM_DIR}\
     --outSAMtype BAM SortedByCoordinate;  
```

This creates the output of a single BAM file where each read is annotated with the header associated with each read in the FASTQ files. The sequencing library is annotated through the **`RG:Z`** tag.

#### Add Cell barcode and UMI to read tags in bam
The cell barcodes and UMIs are already stored in the header of each read. We can use a simple script to loop through each read and move these annotations to the formally-stored tags.

```BASH
(scdiffeq) mvinyard@scdiffeq:~/scripts$ python tag_bam_cb_UMI.py 
```
>Adding CB and UB tags to BAM file...
Number of reads processed: 10000000
2.15 minutes elapsed.
Number of reads barcoded: 10000000 Percentage barcoded: 100.00 %
Number of reads processed: 20000000
4.22 minutes elapsed.
Number of reads barcoded: 20000000 Percentage barcoded: 100.00 %
Number of reads processed: 30000000
6.29 minutes elapsed.
Number of reads barcoded: 30000000 Percentage barcoded: 100.00 %
Number of reads processed: 40000000
8.44 minutes elapsed.
Number of reads barcoded: 40000000 Percentage barcoded: 100.00 %
...
265.17 minutes elapsed.
Number of reads barcoded: 1210000000 Percentage barcoded: 100.00 %
Number of reads processed: 1220000000
267.36 minutes elapsed.
Number of reads barcoded: 1220000000 Percentage barcoded: 100.00 %
Number of reads processed: 1230000000
269.57 minutes elapsed.
Number of reads barcoded: 1230000000 Percentage barcoded: 100.00 %

We can use `samtools view` to peek at the reads. Notice the **`RG`**, **`CB`**, and **`UB`** tags:

```BASH
samtools view LARRY.sorted.tagged.bam | head -2
```
>AACGCCAA-TTGGCGTT:CAGTTA:NS500422_606_HJHGKBGX5_3_21509_14349_14658     16      chr1    3000093 3       2S36M   *       0     0GTTGGCGTTCATTTATTTTTTTTTTTTTTTTTTGGGTG  6//AE/A/A/E/AE/66EEEEEEEEEEEEEEEEAAAAA  NH:i:2  HI:i:1  AS:i:25 nM:i:5  RG:Z:d6_1_2   CB:Z:AACGCCAA-TTGGCGTT   UB:Z:CAGTTA
>
>CTTAGGCC-ATCCCACG:CAAAAA:NB502063_104_HJGTYBGX5_4_22602_17412_17407     256     chr1    3000096 1       4S35M   *       0     0TTAAGCGTATTTCTTTTTTTTTTTTTTTTTTTGGGTGGG AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEA/E/AEE NH:i:4  HI:i:1  AS:i:28 nM:i:3  RG:Z:LK_d6_2_2CB:Z:CTTAGGCC-ATCCCACG   UB:Z:CAAAAA

Next, we index the file:

```BASH
samtools index LARRY.sorted.tagged.bam
```
>LARRY.sorted.tagged.bam.bai


Now we're ready to move to some analysis tools that require bam files as input (i.e., velocyto and Mutect2). We could create bam files for each cell, however 

Create individual bam files for each cell

Index the bam and split into single-cell bam files. 

## Run velocyto
Using the [velocyto CLI](https://velocyto.org/velocyto.py/tutorial/cli.html), I 

```
BAM=/home/mvinyard/data/raw/weinreb_2020/bams/SRR10510898.bam
GTF=/home/mvinyard/ref/mm10/refdata-gex-mm10-2020-A/genes/genes.gtf

velocyto run ${BAM} \
             ${GTF}
```
## Run Mutect2

CRISPR-induced somatic mutations have been used to annotate sub-clonal lineages within the those lineages defined by a static barcode. To determine whether mutational rates in the hematopoeitic system might be sufficient to construct a phylogenetic tree within lineages of the LARRY dataset, we employed *Mutect2*, a sensitive mutational caller to identify mutational profiles in each cell. 

In addition to a static barcode, I wondered if it may be possible to 

## Integration with other processings

**Author pre-processed data**

**CoSpar**

**PRESCIENT**

### [CytoTRACE](https://cytotrace.stanford.edu/)

Using the preprocessed counts matrix available from the Klein Lab Github repository, I ran CytoTRACE following the [CLI tutorial](https://cytotrace.stanford.edu/). I then parsed the CytoTRACE results and incorporated them into an AnnData object. This repo includes a notebook detailing the commands used to run CytoTRACE as well as the subsequent parsing and integration with AnnData.

## Downloadable datasets and useful data subsets

**[Complete LARRY dataset](https://)**

**[LARRY Mutational Profile](https://)**
