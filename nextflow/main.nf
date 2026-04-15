#!/usr/bin/env nextflow
/*
================================================================================
NJU-seq Nextflow Pipeline
================================================================================
NJU-seq: Nm-site Judged Universally sequencing
A comprehensive pipeline for 2'-O-methylation site identification
from RNA-seq data.

Author: NJU-seq Team
Version: 1.0.0
Last Updated: 2026-04-13

Usage:
    nextflow run main.nf [options]

Modules:
    - Bacterial RNA filtering
    - rRNA alignment and analysis
    - mRNA alignment and deduplication
    - Nm-site scoring
    - Visualization and reporting
================================================================================
*/

nextflow.enable.dsl = 2

/*
================================================================================
Include External Modules
================================================================================
*/

include { BACTERIAL_ALIGNMENT } from './modules/bacterial_alignment.nf'
include { BACTERIAL_FILTER } from './modules/bacterial_filter.nf'
include { REMOVE_BACTERIAL_READS } from './modules/remove_bacterial_reads.nf'

include { RRNA_ALIGNMENT } from './modules/rrna_alignment.nf'
include { RRNA_FILTER } from './modules/rrna_filter.nf'
include { RRNA_COUNT } from './modules/rrna_count.nf'
include { RRNA_SCORING } from './modules/rrna_scoring.nf'
include { REMOVE_RRNA_READS } from './modules/remove_rrna_reads.nf'

include { MRNA_ALIGNMENT } from './modules/mrna_alignment.nf'
include { MRNA_FILTER } from './modules/mrna_filter.nf'
include { MRNA_DEDUP } from './modules/mrna_dedup.nf'
include { MRNA_COUNT } from './modules/mrna_count.nf'
include { MRNA_MERGE } from './modules/mrna_merge.nf'

include { MRNA_SCORING } from './modules/mrna_scoring.nf'
include { REGIONAL_ANALYSIS } from './modules/regional_analysis.nf'
include { MOTIF_ANALYSIS } from './modules/motif_analysis.nf'
include { VISUALIZATION } from './modules/visualization.nf'

include { GENERATE_REPORT } from './modules/report.nf'

/*
================================================================================
Main Workflow
================================================================================
*/

workflow NJU_SEQ_PIPELINE {

    take:
    samples_ch
    tissues_ch
    rna_types_ch
    index_ch
    njuseq_dir_ch

    main:

    // Create sample-specific channels
    samples = samples_ch.flatten()

    // =========================================================================
    // Stage 1: Bacterial RNA Filtering
    // =========================================================================

    log.info "Stage 1: Bacterial RNA Filtering"

    bacterial_out = BACTERIAL_ALIGNMENT(
        samples,
        index_ch.bacterial,
        params.data_dir,
        params.output_dir
    )

    bacterial_filter_out = BACTERIAL_FILTER(
        bacterial_out.sam,
        njuseq_dir_ch,
        params.output_dir
    )

    bacterial_remove_out = REMOVE_BACTERIAL_READS(
        bacterial_filter_out.filtered_list,
        samples,
        params.data_dir,
        params.output_dir
    )

    // =========================================================================
    // Stage 2: rRNA Analysis
    // =========================================================================

    log.info "Stage 2: rRNA Alignment and Analysis"

    rrna_out = RRNA_ALIGNMENT(
        bacterial_remove_out.filtered_reads
            .map { it -> tuple(it[0], it[1], it[2]) },
        index_ch.rrna,
        params.output_dir
    )

    rrna_filter_out = RRNA_FILTER(
        rrna_out.sam,
        njuseq_dir_ch,
        params.output_dir
    )

    rrna_count_out = RRNA_COUNT(
        rrna_filter_out.filtered_list,
        samples,
        index_ch.rrna_fasta,
        params.output_dir
    )

    rRNA_by_tissue = rrna_count_out.count
        .map { it -> tuple(it[0], it[1]) }
        .groupTuple()
        .map { it -> tuple(it[0], it[1]) }

    rrna_score_out = RRNA_SCORING(
        rRNA_by_tissue,
        tissues_ch,
        njuseq_dir_ch,
        params.output_dir
    )

    rrna_remove_out = REMOVE_RRNA_READS(
        rrna_filter_out.filtered_list,
        samples,
        params.data_dir,
        params.output_dir
    )

    // =========================================================================
    // Stage 3: mRNA Analysis
    // =========================================================================

    log.info "Stage 3: mRNA Alignment and Analysis"

    mrna_out = MRNA_ALIGNMENT(
        rrna_remove_out.filtered_reads
            .map { it -> tuple(it[0], it[1], it[2]) },
        index_ch.mrna,
        params.output_dir
    )

    mrna_filter_out = MRNA_FILTER(
        mrna_out.sam,
        njuseq_dir_ch,
        params.output_dir
    )

    mrna_dedup_out = MRNA_DEDUP(
        mrna_filter_out.filtered_sam,
        samples,
        index_ch.exon_info,
        params.output_dir
    )

    mrna_count_out = MRNA_COUNT(
        mrna_dedup_out.dedup_sam,
        samples,
        params.data_dir,
        params.output_dir
    )

    mrna_merge_out = MRNA_MERGE(
        mrna_count_out.count,
        samples,
        index_ch.exon_info,
        params.output_dir
    )

    // =========================================================================
    // Stage 4: Scoring and Visualization
    // =========================================================================

    log.info "Stage 4: Scoring and Visualization"

    mrna_by_tissue = mrna_merge_out.merged
        .map { it -> tuple(it[0], it[1]) }
        .groupTuple()
        .map { it -> tuple(it[0], it[1]) }

    mrna_score_out = MRNA_SCORING(
        mrna_by_tissue,
        tissues_ch,
        njuseq_dir_ch,
        params.output_dir
    )

    regional_out = REGIONAL_ANALYSIS(
        mrna_score_out.scored,
        index_ch.transcript_region,
        tissues_ch,
        njuseq_dir_ch,
        params.output_dir
    )

    motif_out = MOTIF_ANALYSIS(
        mrna_score_out.scored,
        index_ch.motif_info,
        tissues_ch,
        njuseq_dir_ch,
        params.output_dir
    )

    viz_out = VISUALIZATION(
        mrna_score_out.scored,
        regional_out.region,
        tissues_ch,
        params.output_dir
    )

    // =========================================================================
    // Stage 5: Report Generation
    // =========================================================================

    log.info "Stage 5: Generating Reports"

    report_out = GENERATE_REPORT(
        mrna_score_out.scored.collect(),
        regional_out.region.collect(),
        motif_out.motif.collect(),
        viz_out.plots.collect(),
        params.output_dir
    )

    emit:
    bacterial_sam = bacterial_out.sam
    rrna_sam = rrna_out.sam
    mrna_sam = mrna_out.sam
    mrna_scored = mrna_score_out.scored
    reports = report_out.report
}

/*
================================================================================
Workflow Execution
================================================================================
*/

workflow {
    // Validate index files exist
    if (!file(params.bac_index).exists()) {
        exit 1, "ERROR: Bacterial index not found: ${params.bac_index}"
    }
    
    if (!file(params.rrna_index).exists()) {
        exit 1, "ERROR: rRNA index not found: ${params.rrna_index}"
    }
    
    if (!file(params.mrna_index).exists()) {
        exit 1, "ERROR: mRNA index not found: ${params.mrna_index}"
    }
    
    if (!file(params.exon_info).exists()) {
        exit 1, "ERROR: Exon info file not found: ${params.exon_info}"
    }
    
    if (!file("${params.index_dir}/ath/transcript_region.info").exists()) {
        exit 1, "ERROR: Transcript region info not found: ${params.index_dir}/ath/transcript_region.info"
    }
    
    if (!file("${params.index_dir}/ath/motif_info.tsv").exists()) {
        exit 1, "ERROR: Motif info file not found: ${params.index_dir}/ath/motif_info.tsv"
    }
    
    // Check if all rRNA fasta files exist
    def rrna_fasta_files = [
        "${params.index_dir}/ath/25s.fa",
        "${params.index_dir}/ath/18s.fa",
        "${params.index_dir}/ath/5-8s.fa"
    ]
    
    def missing_files = []
    for (file_path in rrna_fasta_files) {
        if (!file(file_path).exists()) {
            missing_files << file_path
        }
    }
    
    if (missing_files.size() > 0) {
        exit 1, "ERROR: Missing rRNA fasta files: ${missing_files.join(', ')}"
    }

    // Create input channels
    samples_ch = Channel.fromList(params.samples)
    tissues_ch = Channel.fromList(params.tissues)
    rna_types_ch = Channel.fromList(params.rna_types)
    njuseq_dir_ch = Channel.of(params.njuseq_dir)

    // Create index channel
    index_ch = Channel.from([
        bacterial: params.bac_index,
        rrna: params.rrna_index,
        rrna_fasta: "${params.index_dir}/ath/25s.fa,${params.index_dir}/ath/18s.fa,${params.index_dir}/ath/5-8s.fa",
        mrna: params.mrna_index,
        exon_info: params.exon_info,
        transcript_region: "${params.index_dir}/ath/transcript_region.info",
        motif_info: "${params.index_dir}/ath/motif_info.tsv"
    ])

    // Print workflow information
    log.info """
    ================================================================================
    NJU-seq Pipeline Started
    ================================================================================
    Samples     : ${params.samples.size()}
    Tissues     : ${params.tissues.join(', ')}
    rRNA Types  : ${params.rna_types.join(', ')}
    Output Dir  : ${params.output_dir}
    ================================================================================
    """

    // Execute main pipeline
    NJU_SEQ_PIPELINE(
        samples_ch,
        tissues_ch,
        rna_types_ch,
        index_ch,
        njuseq_dir_ch
    )
}

/*
================================================================================
er Functions
================================================================================
*/

def create_sample_outputs(sample, type, output_dir) {
    return [
        sample: sample,
        type: type,
        output_dir: "${output_dir}/${sample}"
    ]
}

def get_sample_files(sample, data_dir) {
    return [
        r1: "${data_dir}/${sample}/R1.origin.fq.gz",
        r2: "${data_dir}/${sample}/R2.origin.fq.gz"
    ]
}
