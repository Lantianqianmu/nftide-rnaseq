#!/home/zeemeeuw/miniconda3/envs/joint/bin/nextflow


// nextflow rnaseq_pe.nf -with-report nf_rna_report.html -with-timeline nf_rna_timeline.html
// mandatory field: genome, genomeBaseDir, input_csv, -output-dir

params.genome = "hg38"
params.genomeBaseDir = "/home/zeemeeuw/data/ref"
params.gtf = "${params.genomeBaseDir}/${params.genome}/gencode.v39.primary_assembly.annotation.gtf"
params.genomeDir = "${params.genomeBaseDir}/${params.genome}/${params.genome}-STAR/"
params.input_csv = '/home/zeemeeuw/data/nextflow/data/rnaseq/DS601_RNA/samplesheet.csv'


// Channel
//     .fromFilePairs(params.reads, checkIfExists: true, flat: true)
//     //.view()
//     .set { read_pairs_ch }

// Create input channel from the contents of a CSV file

process CUTADAPT {
    tag "cutadapt on ${id}"
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
    tag "STAR on ${id}"
    // publishDir "${params.workingDir}/${id}/STAR", mode: 'copy', overwrite: false

    input:
    tuple val(id), path(read1), path(read2)

    output:
    val(id), emit: aligned_sample
    tuple val(id), path("*_aligned.bam"), emit: aligned_bam
    tuple val(id), path("*_aligned_filtered.bam"), emit: aligned_filtered_bam
    path("*_aligned_filtered.bam"), emit: bams
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
    tag "featureCounts on all samples"
    // publishDir "${params.workingDir}", mode: 'copy', overwrite: true

    input:
    path bams

    output:
    path("count_matrix.txt"), emit: count_matrix
    path("featureCounts.log"), emit: featureCounts_log

    script:
    """
    featureCounts \
    -p \
    --primary \
    -T ${task.cpus} \
    -t gene -g gene_id \
    -a ${params.gtf} -o count_matrix.txt \
    ${bams} > "featureCounts.log"
    """

}



workflow {
    main:

    read_pairs_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)] }

    log.info """\
      nftide-rnaseq
      ===================================
      genome     : ${params.genome}
      gtf        : ${params.gtf}
      genomeDir  : ${params.genomeDir}
      projectDir : ${projectDir}
      workingDir : ${workflow.outputDir}
    """.stripIndent()

    CUTADAPT(read_pairs_ch)
    // CUTADAPT.out.trimmed_reads.view()
    STAR(CUTADAPT.out.trimmed_reads)
    FEATURECOUNTS(STAR.out.bams.collect())

    publish:
    cutadapt_fastqs = CUTADAPT.out.trimmed_reads
    cutadapt_logs = CUTADAPT.out.cutadapt_log
    star_aligned_bams = STAR.out.aligned_bam
    star_aligned_filtered_bams = STAR.out.aligned_filtered_bam
    star_qc = STAR.out.aligned_log
    featureCounts_logs = FEATURECOUNTS.out.featureCounts_log
    featureCounts_mat = FEATURECOUNTS.out.count_matrix

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
}


