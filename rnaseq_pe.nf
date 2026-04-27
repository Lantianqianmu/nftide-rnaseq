#!/home/zeemeeuw/miniconda3/envs/joint/bin/nextflow


// nextflow rnaseq_pe.nf -with-report nf_rna_report.html -with-timeline nf_rna_timeline.html
// mandatory field: genome, genomeBaseDir, input_csv, -output-dir

// nextflow run rnaseq_pe.nf -output-dir /data/xrz/charseq/rnaseq_nf_out -with-report nf_rna_report.html -with-timeline nf_rna_timeline.html 
params.gtf = "/data/xrz/ref/hg38/gencode.v39.primary_assembly.annotation.gtf"
params.genomeDir = "/data/xrz/ref/hg38/hg38-STAR/"
params.input_csv = '/data/xrz/charseq/nftide-rnaseq/samplesheet.csv'


// Channel
//     .fromFilePairs(params.reads, checkIfExists: true, flat: true)
//     //.view()
//     .set { read_pairs_ch }

// Create input channel from the contents of a CSV file


process CUTADAPT {
    tag "cutadapt on ${id}..."
    // publishDir "${params.workingDir}/${id}/cutadapt", mode: 'copy', overwrite: false

    input:
    tuple val(id), path(read1), path(read2)

    output:
    tuple val(id), path("*_cutadapt_R1.fq.gz"), path("*_cutadapt_R2.fq.gz"), emit: trimmed_reads
    tuple val(id), path("*_cutadapt.log"), emit: cutadapt_log

    script:
    """
    cutadapt \
        -j ${task.cpus} -m 15:15 \
        -a "CTGTCTCTTATACACATCT" \
        -A "CTGTCTCTTATACACATCT" \
        --pair-filter=any \
        -o ${id}_cutadapt_R1.fq.gz \
        -p ${id}_cutadapt_R2.fq.gz \
        ${read1} \
        ${read2} > "${id}_cutadapt.log"

    """
}

process STAR {
    tag "STAR on ${id}..."
    // publishDir "${params.workingDir}/${id}/STAR", mode: 'copy', overwrite: false

    input:
    tuple val(id), path(read1), path(read2)

    output:
    val(id), emit: aligned_sample
    tuple val(id), path("*_aligned.bam"), emit: aligned_bam
    tuple val(id), path("*_aligned_filtered.bam"), emit: aligned_filtered_bam
    // path("*_aligned_filtered.bam"), emit: bams
    tuple val(id), path("*_Log.final.out"), emit: aligned_log

    script:
    """
    STAR \
        --runMode alignReads \
        --runThreadN ${task.cpus} \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNoverLmax 0.1 \
        --limitOutSJcollapsed 2000000 \
        --chimOutType WithinBAM \
        --readFilesIn ${read1} ${read2} \
        --genomeDir ${params.genomeDir} \
        --outFileNamePrefix ${id}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN ${task.cpus} 

    mv ${id}_Aligned.sortedByCoord.out.bam ${id}_aligned.bam

    samtools view -@ ${task.cpus} -h -b -F 772 -q 30 ${id}_aligned.bam > ${id}_aligned_filtered.bam
    """

}

process FEATURECOUNTS {
    tag "featureCounts on all samples..."
    // publishDir "${params.workingDir}", mode: 'copy', overwrite: true

    input:
    val order
    val bams

    output:
    path("featureCounts_matrix.txt"), emit: count_matrix
    path("featureCounts.log"), emit: featureCounts_log

    script:
    def sorted_bams = bams.sort { bamTuple ->
        order.indexOf(bamTuple[0])
    }
    
    def bam_paths = sorted_bams.collect { item -> item[1] }
    

    def clean_paths = bam_paths
        .collect { item -> item.toString() }
        .collect { item -> item.replaceAll("[\\[\\]',]", "").trim() }
        .findAll { item -> item }
        .join(' ')

    """
    BAM_ARRAY=()
    for bam in ${clean_paths}; do

        base_name=\$(basename "\$bam")
        BAM_ARRAY+=("\$base_name")
        ln -sf "\$bam" "\$base_name"
    done

    featureCounts \
    -p \
    --primary \
    -T ${task.cpus} \
    -t gene -g gene_id \
    -a ${params.gtf} -o featureCounts_matrix.txt \
    "\${BAM_ARRAY[@]}" 2> "featureCounts.log"

    """
}

process MAKEMATRIX {
    tag "Making the final count matrix..."
    // publishDir "${params.workingDir}", mode: 'copy', overwrite: true

    input:
    val count_mat

    output:
    path("count_matrix.txt"), emit: final_count_matrix

    script:

    """
        cat ${count_mat} | perl -alne '
            BEGIN{
                \$, = "\t";
                open IN, "<", shift;
                while(<IN>){
                    next if (/^#/);
                    if(/gene_id "(.*?)"/){\$id = \$1;}
                    if(/gene_type "(.*?)"/){\$type = \$1;}
                    if(/gene_name "(.*?)"/){\$name = \$1;}
                    \$gname{\$id} = \$name; \$gtype{\$id} = \$type; 
                }
                close IN;
            }
            next if (\$. == 1);
            if (\$. == 2){
                foreach my \$str (@F) {
                    \$str =~ s/_aligned_filtered\\.bam//;
                }
                print "Genename", "Genetype", @F;
            }else{
                print(\$gname{\$F[0]}, \$gtype{\$F[0]}, @F);
            }
        ' ${params.gtf} > count_matrix.txt

    """
}



workflow {
    main:
    read_pairs_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)] }
    sample_order = read_pairs_ch.map { sample, _f1, _f2 -> sample }.collect()
    

    log.info """\
      nftide-rnaseq
      ===================================
      gtf        : ${params.gtf}
      genomeDir  : ${params.genomeDir}
      projectDir : ${projectDir}
      workingDir : ${workflow.outputDir}
    """.stripIndent()

    CUTADAPT(read_pairs_ch)
    STAR(CUTADAPT.out.trimmed_reads)

    FEATURECOUNTS(sample_order, STAR.out.aligned_filtered_bam.collect(flat: false))

    MAKEMATRIX(FEATURECOUNTS.out.count_matrix)

    publish:
    cutadapt_fastqs = CUTADAPT.out.trimmed_reads
    cutadapt_logs = CUTADAPT.out.cutadapt_log
    star_aligned_bams = STAR.out.aligned_bam
    star_aligned_filtered_bams = STAR.out.aligned_filtered_bam
    star_qc = STAR.out.aligned_log
    featureCounts_logs = FEATURECOUNTS.out.featureCounts_log
    featureCounts_mat = FEATURECOUNTS.out.count_matrix
    final_mat = MAKEMATRIX.out.final_count_matrix

}


output {
    cutadapt_fastqs {
        path { sample, _f1, _f2 -> "${sample}/cutadapt" }
    }
    cutadapt_logs {
        path { sample, _f1 -> "${sample}/cutadapt" }
    }
    star_aligned_bams {
        path { sample, _f1 -> "${sample}/STAR" }
    }
    star_aligned_filtered_bams {
        path { sample, _f1 -> "${sample}/STAR" }
    }
    star_qc {
        path { sample, _f1 -> "${sample}/STAR" }
    }
    featureCounts_logs {
    }    
    featureCounts_mat {
    }
    final_mat {
    }
}


