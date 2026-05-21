#!/usr/bin/env nextflow

/*
 * ============================================================
 *  TERRA-QUANT: Telomeric Repeat-containing RNA Quantification
 * ============================================================
 *  Quantifies TERRA expression from bulk RNA-seq using the
 *  T2T-CHM13 reference with custom TERRA/ITS annotations.
 *
 *  Usage:
 *    nextflow run main.nf --samplesheet samples.csv --star_index /path/to/star_index
 *
 *    # Build STAR index on-the-fly:
 *    nextflow run main.nf --samplesheet samples.csv --genome_fasta /path/to/chm13.fa
 *
 *  Samplesheet CSV format (no header):
 *    sample_id,fastq_r1,fastq_r2
 *
 *  Multiple rows per sample_id are allowed (lanes/runs will be merged).
 */

nextflow.enable.dsl = 2

// ---- Parameter validation ----
if (!params.samplesheet) {
    error "ERROR: --samplesheet is required. Provide a CSV with columns: sample,fastq_1,fastq_2"
}
if (!params.star_index && !params.genome_fasta) {
    error "ERROR: Provide either --star_index (pre-built) or --genome_fasta (to build index)."
}

// ---- Log ----
log.info """
=============================================================
 TERRA-QUANT  v1.0
=============================================================
 samplesheet  : ${params.samplesheet}
 star_index   : ${params.star_index ?: 'will be built from genome_fasta'}
 gtf          : ${params.gtf}
 strandedness : ${params.strandedness}
 outdir       : ${params.outdir}
 run_bbduk    : ${params.run_bbduk}
 run_yarn     : ${params.run_yarn}
=============================================================
"""

// ============================================================
//  PROCESSES
// ============================================================

process GUNZIP_GTF {
    label 'low'

    input:
    path gtf_gz

    output:
    path "*.gtf", emit: gtf

    script:
    """
    gunzip -c ${gtf_gz} > ${gtf_gz.baseName}
    """
}

process STAR_INDEX {
    label 'high_mem'
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    """
    mkdir star_index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang 100
    """
}

process TRIM_GALORE {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed/${sample_id}", mode: 'copy', pattern: '*_trimming_report.txt'
    publishDir "${params.outdir}/trimmed/${sample_id}", mode: 'copy', pattern: '*_fastqc.*'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: reports
    path "*_fastqc.*",            emit: fastqc

    script:
    """
    trim_galore \\
        --paired \\
        --quality ${params.trim_quality} \\
        --fastqc \\
        --illumina \\
        --gzip \\
        --cores ${task.cpus} \\
        ${r1} ${r2}
    """
}

process CAT_FASTQ {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(r1_files), path(r2_files)

    output:
    tuple val(sample_id), path("${sample_id}_merged_R1.fq.gz"), path("${sample_id}_merged_R2.fq.gz"), emit: merged_reads

    script:
    """
    cat ${r1_files} > ${sample_id}_merged_R1.fq.gz
    cat ${r2_files} > ${sample_id}_merged_R2.fq.gz
    """
}

process STAR_ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: '*.bam'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path star_index

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample_id}.Log.final.out", emit: log

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --runMode alignReads \\
        --readFilesCommand zcat \\
        --readFilesIn ${r1} ${r2} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${sample_id}.
    """
}

process SAMTOOLS_INDEX {
    tag "${sample_id}"
    publishDir "${params.outdir}/alignments", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: bam

    script:
    """
    samtools index ${bam}
    """
}

process SAMTOOLS_NAME_SORT {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.namesorted.bam"), emit: bam

    script:
    """
    samtools sort -n -@ ${task.cpus} ${bam} -o ${sample_id}.namesorted.bam
    """
}

process HTSEQ_COUNT {
    tag "${sample_id}"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path gtf

    output:
    path "${sample_id}.count.txt", emit: counts

    script:
    """
    htseq-count \\
        -f bam \\
        -r name \\
        -s ${params.strandedness} \\
        -t exon \\
        --idattr gene_name \\
        -m intersection-nonempty \\
        --nonunique all \\
        ${bam} \\
        ${gtf} > ${sample_id}.count.txt
    """
}

process BAMCOVERAGE {
    tag "${sample_id}"
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "*.bigwig", emit: bigwigs

    script:
    """
    bamCoverage \\
        -b ${bam} \\
        --normalizeUsing RPKM \\
        -p ${task.cpus} \\
        -o ${sample_id}.RPKM.unstranded.bigwig \\
        -of bigwig

    bamCoverage \\
        -b ${bam} \\
        --filterRNAstrand forward \\
        --normalizeUsing RPKM \\
        -p ${task.cpus} \\
        -o ${sample_id}.RPKM.forward.bigwig \\
        -of bigwig

    bamCoverage \\
        -b ${bam} \\
        --filterRNAstrand reverse \\
        --normalizeUsing RPKM \\
        -p ${task.cpus} \\
        -o ${sample_id}.RPKM.reverse.bigwig \\
        -of bigwig
    """
}

process BBDUK_TELO {
    tag "${sample_id}"
    publishDir "${params.outdir}/bbduk_telo", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path telo_ref

    output:
    path "${sample_id}_telo_content.stats.txt", emit: stats
    path "${sample_id}_telo_content*.fa",       emit: telo_reads

    script:
    """
    bbduk.sh \\
        overwrite=t \\
        in=${r1} \\
        in2=${r2} \\
        ref=${telo_ref} \\
        k=24 hdist=2 \\
        threads=${task.cpus} \\
        outm=${sample_id}_telo_content_R1.fa \\
        outm2=${sample_id}_telo_content_R2.fa \\
        stats=${sample_id}_telo_content.stats.txt
    """
}

process SUMMARIZE_COUNTS {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path count_files

    output:
    path "raw_counts_matrix.csv",        emit: raw_counts
    path "TERRA_repeat_counts.csv",      emit: terra_repeat
    path "TERRA_subtelo_counts.csv",     emit: terra_subtelo
    path "gene_counts_no_TERRA.csv",     emit: gene_counts

    script:
    """
    summarize_counts.R . .
    """
}

process YARN_NORMALIZE {
    publishDir "${params.outdir}/normalized", mode: 'copy'

    input:
    path gene_counts
    path terra_counts
    path annotation

    output:
    path "YARN_normalized_count_withTERRA.csv", emit: norm_all
    path "YARN_normalized_TERRA_count.csv",     emit: norm_terra
    path "*.png",                               emit: plots

    script:
    """
    yarn_normalize.R \\
        --counts ${gene_counts} \\
        --terra_counts ${terra_counts} \\
        --annotation ${annotation} \\
        --target ${params.yarn_target} \\
        --outdir .
    """
}

// ============================================================
//  WORKFLOW
// ============================================================

workflow {

    // ---- Prepare GTF ----
    gtf_file = file(params.gtf, checkIfExists: true)
    if (params.gtf.endsWith('.gz')) {
        GUNZIP_GTF(gtf_file)
        ch_gtf = GUNZIP_GTF.out.gtf
    } else {
        ch_gtf = Channel.value(gtf_file)
    }

    // ---- Prepare STAR index ----
    if (params.star_index) {
        ch_star_index = Channel.value(file(params.star_index, checkIfExists: true))
    } else {
        ch_genome = Channel.value(file(params.genome_fasta, checkIfExists: true))
        STAR_INDEX(ch_genome, ch_gtf)
        ch_star_index = STAR_INDEX.out.index
    }

    // ---- Parse samplesheet ----
    // Format: sample,fastq_1,fastq_2
    ch_reads = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample = row.sample
            def r1 = file(row.fastq_1, checkIfExists: true)
            def r2 = file(row.fastq_2, checkIfExists: true)
            tuple(sample, r1, r2)
        }

    // ---- Step 1: Trim (per FASTQ pair) ----
    TRIM_GALORE(ch_reads)

    // ---- Step 2: Merge trimmed FASTQs per sample (if multiple pairs) ----
    ch_trimmed_grouped = TRIM_GALORE.out.trimmed_reads
        .groupTuple()
        .map { sample_id, r1_list, r2_list ->
            tuple(sample_id, r1_list.flatten(), r2_list.flatten())
        }

    // Only merge if needed, but CAT_FASTQ handles single pairs fine too
    CAT_FASTQ(ch_trimmed_grouped)

    // ---- Step 3: Align ----
    STAR_ALIGN(CAT_FASTQ.out.merged_reads, ch_star_index)

    // ---- Step 3b: Index BAM ----
    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)

    // ---- Step 3c: Name-sort BAM for htseq-count (avoids mate-pair buffer memory) ----
    SAMTOOLS_NAME_SORT(STAR_ALIGN.out.bam)

    // ---- Step 4: Count ----
    HTSEQ_COUNT(SAMTOOLS_NAME_SORT.out.bam, ch_gtf)

    // ---- Step 5: Summarize counts ----
    SUMMARIZE_COUNTS(HTSEQ_COUNT.out.counts.collect())

    // ---- Optional: bamCoverage ----
    if (params.run_bamcoverage) {
        BAMCOVERAGE(SAMTOOLS_INDEX.out.bam)
    }

    // ---- Optional: BBDuk telomeric content ----
    if (params.run_bbduk) {
        ch_telo_ref = Channel.value(file(params.telo_ref, checkIfExists: true))
        BBDUK_TELO(CAT_FASTQ.out.merged_reads, ch_telo_ref)
    }

    // ---- Optional: YARN normalization ----
    if (params.run_yarn) {
        if (!params.annotation || !params.yarn_target) {
            error "ERROR: --annotation and --yarn_target required when --run_yarn is true"
        }
        ch_annotation = Channel.value(file(params.annotation, checkIfExists: true))
        YARN_NORMALIZE(
            SUMMARIZE_COUNTS.out.gene_counts,
            SUMMARIZE_COUNTS.out.terra_repeat,
            ch_annotation
        )
    }
}
