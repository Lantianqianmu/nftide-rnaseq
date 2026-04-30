# nftide-rnaseq #
Nextflow pipeline for pair-end RNA-seq.

## Introduction ##
The pipeline uses cutadapt to remove nextera adaptors, aligns fastqs with STAR and quantifies with featureCounts.  

## Software dependencies ##
Dependencies  | Version
------------- | -------------
nextflow | 25.10
openjdk | 17.0.8
Perl | 5.32.1
python | 3.10.12
cutadapt | 4.5
samtools | 1.18
star | 2.7.11a
subread | 2.0.6

## Installation ##
(1) Create a conda environment with 
```
mamba create -n nftide-rnaseq nextflow=25.10.2 python cutadapt samtools star subread
``` 
and activate the environment with  
```
mamba activate nftide-rnaseq
```

(2) Clone the repository with `git clone`, and execute
```
cd nftide-rnaseq
```
(3) You also need to prepare STAR index files for your genome. Please refer to their manuals to generate indexed genome files.

## Usage ##
(1) Prepare the `samplesheet.csv`. The csv file __must__ contain 3 columns with defined column names:  
`sample`: Name of the sequenced library. For example, `demo-1`. It will be the prefix of the output. Note: Different fastqs with same sample name will be merged before processing.  
`fastq_1`: Path to read 1.  
`fastq_2`: Path to read 2.  

(2) Run nextflow pipeline.
```
nextflow run rnaseq_pe.nf \
  -output-dir outdir \
  --gtf gtffile \
  --genomeDir STARindexfolder \
  --input_csv samplesheet.csv \
  -with-report outdir/nf_rna_report.html \
  -with-timeline outdir/nf_rna_timeline.html \
  -bg
```
`-output-dir`: Path to the output directory.  
`--input_csv`: Path to samplesheet.csv as described in **step (1)**.  
`--gtf`: gtf annotation file for STAR.  
`--genomeDir`: STAR index folder.  
By default, the pipeline allows 2 samples to be processed in parallel. To change this behavior, modify _maxForks_ in __nextflow.config__.

## Expected output ##
The pipeline creates subfolders (named by samples in the samplesheet) in -output-dir. In each subfolder, there will be a __cutadapt__ and a __STAR__ folder.  
The count matrix is outputed as __count_matrix.txt__ in -output-dir. The column order of __count_matrix.txt__ is determined by the "sample" column of input csv.
