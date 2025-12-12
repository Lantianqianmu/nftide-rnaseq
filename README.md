# nftide-rnaseq #
Nextflow pipeline for pair-end RNA-seq.

## Introduction ##
The pipeline uses cutadapt to remove nextera adaptors, aligns fastqs with STAR and quantifies with featureCounts.  

## Software dependencies ##
Dependencies  | Version
------------- | -------------
nextflow | 25.10
Perl | 5.32.1
python | 3.10.12
cutadapt | 4.5
samtools | 1.18
star | 2.7.11a
subread | 2.0.6

## Usage ##
(1) Create a conda environment with 
```
conda create -n nftide-rnaseq nextflow=25.10.2 python cutadapt samtools star subread
``` 
and activate the environment with  
```
conda activate nftide-rnaseq
```

(2) Clone the repository with `git clone`, and execute
```
cd nftide-rnaseq
```

(3) Run nextflow pipeline.
```
nextflow rnaseq_pe.nf \
  -output-dir outdir \
  --gtf gtffile \
  --genomeDir STARindexfolder \
  --input_csv samplesheet.csv \
  -with-report outdir/nf_rna_report.html \
  -with-timeline outdir/nf_rna_timeline.html
```
where __samplesheet.csv__ must have 3 columns named "sample", "fastq_1" and "fastq_2".  
By default, the pipeline allows 2 samples to be processed in parallel. To change this behavior, modify maxForks in __nextflow.config__.

## Expected output ##
The pipeline creates subfolders (named by samples in the samplesheet) in -output-dir. In each subfolder, there will be a __cutadapt__ and a __STAR__ folder.  
The count matrix is outputed as __count_matrix.txt__ in -output-dir. The column order of __count_matrix.txt__ is determined by the "sample" column of input csv.
