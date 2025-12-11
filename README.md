# nftide-rnaseq #
Nextflow pipeline for RNA-seq.

### Software dependencies ###
Dependencies  | Version
------------- | -------------
Perl | 5.32.1
python | 3.10.12
cutadapt | 4.5
pigz | 2.6
samtools | 1.18
star | 2.7.11a

### Usage ###

```
nextflow rnaseq_pe.nf \
  --gtf gtffile \
  --genomeDir STARindexfolder \
  --input_csv samplesheet.csv \
  -with-report nf_rna_report.html \
  -with-timeline nf_rna_timeline.html
```

where __samplesheet.csv__ must have 3 columns named "sample", "fastq_1" and "fastq_2".



