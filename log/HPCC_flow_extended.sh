#!/usr/bin/env bash
################################################################################
# HPCC_flow_extended.sh - Extended Arabidopsis thaliana Analysis Pipeline
#
# Description:
#   This extended pipeline includes comprehensive analysis beyond basic Nm-site
#   identification, featuring:
#     - Regional distribution analysis (UTR, CDS, intron)
#     - Alternative splicing site analysis
#     - Gene ontology enrichment
#     - Sequence logo generation
#     - Venn diagram comparison
#     - Statistical analysis
#
# Author: NJU-seq Team
# Last Updated: 2026-04-13
################################################################################

set -euo pipefail

################################################################################
# Configuration Section
################################################################################

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_ROOT="${SCRIPT_DIR}/.."
readonly NJU_SEQ_DIR="${PROJECT_ROOT}"  # 修正：项目根目录就是NJU_seq

readonly SAMPLES=(NJU6220 NJU6221 NJU6222 NJU6223 \
                  NJU6224 NJU6225 NJU6226 NJU6227 \
                  NJU6228 NJU6229 NJU6230 NJU6231 \
                  NJU6232 NJU6233 NJU6234 NJU6235)

readonly RNAS=(25s 18s 5-8s)
readonly READS=(R1 R2)

readonly THREADS=24
readonly BOWTIE_THREADS=22
readonly PARALLEL_JOBS=20

readonly DATA_DIR="${PROJECT_ROOT}/data"  # 修正：使用实际存在的data目录
readonly OUTPUT_DIR="${PROJECT_ROOT}/NJU_seq_output"  # 保持不变
readonly INDEX_DIR="${PROJECT_ROOT}/index"  # 保持不变

readonly ANALYSIS_DIR="${OUTPUT_DIR}/extended_analysis"
readonly FIGURES_DIR="${ANALYSIS_DIR}/figures"
readonly STATS_DIR="${ANALYSIS_DIR}/statistics"
readonly REPORTS_DIR="${ANALYSIS_DIR}/reports"

################################################################################
# Function Definitions
################################################################################

log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $*"
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

check_dependencies() {
    log_info "Checking dependencies..."
    
    # 检查必要的目录
    for dir in "${NJU_SEQ_DIR}" "${OUTPUT_DIR}" "${INDEX_DIR}"; do
        if [[ ! -d "${dir}" ]]; then
            log_error "Directory not found: ${dir}"
            log_error "Please ensure all required directories exist before running the pipeline"
            exit 1
        fi
    done
    
    # 检查必要的文件
    local required_files=(
        "${NJU_SEQ_DIR}/mrna_analysis/main_transcript_3.pl"
        "${NJU_SEQ_DIR}/mrna_analysis/main_transcript_4.pl"
        "${NJU_SEQ_DIR}/mrna_analysis/exon_distance_2.pl"
        "${NJU_SEQ_DIR}/mrna_analysis/motif_nm.pl"
        "${NJU_SEQ_DIR}/presentation/point_distribution.R"
        "${NJU_SEQ_DIR}/presentation/exon_distance.R"
        "${NJU_SEQ_DIR}/presentation/seq_logo.sh"
        "${INDEX_DIR}/ath/transcript_region.info"
        "${INDEX_DIR}/ath/exon_distance.info"
        "${INDEX_DIR}/ath/ath.fa"
    )
    
    local missing_files=()
    for file in "${required_files[@]}"; do
        if [[ ! -f "${file}" ]]; then
            missing_files+=("${file}")
        fi
    done
    
    if [[ ${#missing_files[@]} -gt 0 ]]; then
        log_error "Missing required files:"
        for file in "${missing_files[@]}"; do
            log_error "  - ${file}"
        done
        log_error "Please ensure all required files exist before running the pipeline"
        exit 1
    fi
    
    log_info "All dependencies checked successfully"
}

setup_directories() {
    log_info "Creating extended analysis directory structure..."
    mkdir -p "${ANALYSIS_DIR}"
    mkdir -p "${FIGURES_DIR}"
    mkdir -p "${STATS_DIR}"
    mkdir -p "${REPORTS_DIR}"
    log_info "Directory structure created"
}

################################################################################
# Step 5: Regional Distribution Analysis
################################################################################

step5_region_distribution() {
    log_info "=== Step 5: Regional Distribution Analysis ==="

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        log_info "Analyzing regional distribution for ${tissue}..."

        perl "${NJU_SEQ_DIR}"/mrna_analysis/main_transcript_3.pl \
            "${INDEX_DIR}"/ath/ath_transcript_region.tsv \
            <(sed '1d' "${OUTPUT_DIR}"/mrna.filter."${tissue}"_RF.tsv) \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_region_tmp.tsv |
            perl "${NJU_SEQ_DIR}"/mrna_analysis/main_transcript_4.pl \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_region_norm.tsv &

    done
    wait

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        Rscript "${NJU_SEQ_DIR}"/presentation/point_distribution.R \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_region_norm.tsv \
            "${FIGURES_DIR}"/mrna.filter."${tissue}"_region.pdf \
            --width 48 --height 18 --format pdf --dpi 300
    done

    log_info "Regional distribution analysis completed"
}

step5_filter_by_region() {
    log_info "Filtering Nm sites by genomic region..."

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        # 使用region_norm.tsv文件，因为region.tsv不存在
        perl "${NJU_SEQ_DIR}"/mrna_analysis/filter_transcript.pl \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_region_norm.tsv \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_transcript.tsv &

    done
    wait

    log_info "Region filtering completed"
}

step5_exon_intron_analysis() {
    log_info "Analyzing exon/intron distribution..."

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        log_info "Analyzing exon distance for ${tissue}..."

        perl "${NJU_SEQ_DIR}"/mrna_analysis/exon_distance_2.pl \
            "${INDEX_DIR}"/ath/exon_distance.info \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_region_tmp.tsv \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_exon_site_bar.tsv \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_exon_site_porta.tsv \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_exon_site.tsv &
    done
    wait

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        Rscript "${NJU_SEQ_DIR}"/presentation/exon_distance.R \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_exon_site_bar.tsv \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_exon_site_porta.tsv \
            "${FIGURES_DIR}"/mrna.filter."${tissue}"_exon_site.pdf
    done

    log_info "Exon/intron analysis completed"
}

################################################################################
# Step 6: Alternative Splicing Analysis
################################################################################

step6_identify_alter_splice() {
    log_info "=== Step 6: Alternative Splicing Analysis ==="
    log_info "Identifying alternative splicing sites..."

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        log_info "Identifying alternative splicing sites for ${tissue}..."

        perl "${NJU_SEQ_DIR}"/mrna_analysis/judge_altersplice.pl \
            --transwording "mRNA" \
            --geneid "gene_id=" \
            --alter "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_alter.tsv \
            --unique "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_unique.tsv \
            <index/ath/annotation_processed.gff3 &

    done
    wait

    log_info "Alternative splicing identification completed"
}

step6_judge_splice_exon() {
    log_info "Judging splice exon usage..."

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        perl "${NJU_SEQ_DIR}"/presentation/judge_splice_exon.pl \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_alter.tsv \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_constant.tsv \
            "${OUTPUT_DIR}"/mrna.filter."${tissue}"_RF.tsv \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_splice_junction.tsv &

    done
    wait

    log_info "Splice exon analysis completed"
}

step6_splice_statistics() {
    log_info "Generating splice junction statistics..."

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        perl "${NJU_SEQ_DIR}"/mrna_analysis/stat_altersplice_1.pl \
            "${INDEX_DIR}"/ath/ath_altersplice_gene.info \
            <"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_splice_junction.tsv \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_splice_stats.tsv &

    done
    wait

    perl "${NJU_SEQ_DIR}"/presentation/splice_pie.R \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_leaf_splice_stats.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_root_splice_stats.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_stem_splice_stats.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_flower_splice_stats.tsv \
        "${FIGURES_DIR}"/splice_junction_comparison.pdf

    log_info "Splice statistics completed"
}

################################################################################
# Step 7: Gene Expression and Annotation Analysis
################################################################################

step7_gene_expression() {
    log_info "=== Step 7: Gene Expression and Annotation Analysis ==="

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        log_info "Analyzing gene expression for ${tissue}..."

        perl "${NJU_SEQ_DIR}"/mrna_analysis/add_gene_name.pl \
            --id "gene_id=" \
            --name "gene_name=" \
            --file "${INDEX_DIR}"/ath/ath_gene_info.tsv \
            --col 9 \
            <(sed '1d' "${OUTPUT_DIR}"/mrna.filter."${tissue}"_RF.tsv) \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_with_gene.tsv &

    done
    wait

    log_info "Gene expression analysis completed"
}

step7_transcript_selection() {
    log_info "Selecting representative transcripts..."

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        perl "${NJU_SEQ_DIR}"/mrna_analysis/main_transcript_1.pl \
            --geneid "gene_id=" \
            --transid "transcript_id=" \
            <"${INDEX_DIR}"/ath/ath_transcript_info.gtf \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_representative_transcripts.tsv &

    done
    wait

    log_info "Transcript selection completed"
}

step7_filter_overlapping_genes() {
    log_info "Filtering overlapping gene Nm sites..."

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        perl "${NJU_SEQ_DIR}"/mrna_analysis/filter_overlapgene.pl \
            <"${INDEX_DIR}"/ath/ath_gene_annotation.gff3 \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_nonoverlap.tsv &

    done
    wait

    log_info "Overlapping gene filtering completed"
}

################################################################################
# Step 8: Motif and Sequence Analysis
################################################################################

step8_motif_analysis() {
    log_info "=== Step 8: Motif and Sequence Analysis ==="

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        log_info "Analyzing sequence motifs for ${tissue}..."

        perl "${NJU_SEQ_DIR}"/mrna_analysis/motif_nm.pl \
            "${INDEX_DIR}"/ath/ath.fa \
            <(sed '1d' "${OUTPUT_DIR}"/mrna.filter."${tissue}"_RF.tsv) \
            10 10 \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_motif.tsv &

    done
    wait

    log_info "Motif analysis completed"
}

step8_motif_comparison() {
    log_info "Comparing motifs across tissues..."

    perl "${NJU_SEQ_DIR}"/mrna_analysis/motif_compare.pl \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_leaf_motif.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_root_motif.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_stem_motif.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_flower_motif.tsv \
        >"${ANALYSIS_DIR}"/motif_comparison.tsv

    log_info "Motif comparison completed"
}

step8_sequence_logo() {
    log_info "Generating sequence logos..."

    perl "${NJU_SEQ_DIR}"/presentation/seq_logo.sh \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_leaf_motif.tsv \
        "${INDEX_DIR}"/ath/ath.fa \
        "${FIGURES_DIR}"/sequence_logo_leaf.pdf

    perl "${NJU_SEQ_DIR}"/presentation/seq_logo.sh \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_root_motif.tsv \
        "${INDEX_DIR}"/ath/ath.fa \
        "${FIGURES_DIR}"/sequence_logo_root.pdf

    perl "${NJU_SEQ_DIR}"/presentation/seq_logo.sh \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_stem_motif.tsv \
        "${INDEX_DIR}"/ath/ath.fa \
        "${FIGURES_DIR}"/sequence_logo_stem.pdf

    perl "${NJU_SEQ_DIR}"/presentation/seq_logo.sh \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_flower_motif.tsv \
        "${INDEX_DIR}"/ath/ath.fa \
        "${FIGURES_DIR}"/sequence_logo_flower.pdf

    log_info "Sequence logo generation completed"
}

step8_codon_analysis() {
    log_info "=== Step 8: Codon Context Analysis ==="
    log_info "Analyzing codon context..."

    log_info "Extracting annotation data once..."
    # 直接解压到临时文件，不重新压缩
    gunzip -c "${INDEX_DIR}"/ath/ath.gff3.gz 2>/dev/null >"${ANALYSIS_DIR}"/ath.gff3.tmp &
    gunzip -c "${INDEX_DIR}"/ath/ath.gff3.gz 2>/dev/null >"${ANALYSIS_DIR}"/ath_cds.gff3.tmp &
    wait

    log_info "Generating codon information files..."
    cat "${ANALYSIS_DIR}"/ath_cds.gff3.tmp \
        | awk '$3=="CDS"{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
        | perl "${NJU_SEQ_DIR}"/mrna_analysis/codon_distance_1.pl \
            "${INDEX_DIR}"/ath/ath_main_transcript.txt \
            >"${ANALYSIS_DIR}"/ath_main_transcript_start_codon.tsv &

    cat "${ANALYSIS_DIR}"/ath_cds.gff3.tmp \
        | awk '$3=="CDS"{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
        | perl "${NJU_SEQ_DIR}"/mrna_analysis/codon_distance_1.pl \
            "${INDEX_DIR}"/ath/ath_main_transcript.txt \
            >"${ANALYSIS_DIR}"/ath_main_transcript_stop_codon.tsv &

    cat "${ANALYSIS_DIR}"/ath.gff3.tmp \
        | awk '$3=="transcript"||$3=="exon"' \
        | perl "${NJU_SEQ_DIR}"/mrna_analysis/codon_distance_2.pl \
            "${ANALYSIS_DIR}"/ath_main_transcript_start_codon.tsv \
            >"${ANALYSIS_DIR}"/ath_main_transcript_start_codon.yml &

    cat "${ANALYSIS_DIR}"/ath.gff3.tmp \
        | awk '$3=="transcript"||$3=="exon"' \
        | perl "${NJU_SEQ_DIR}"/mrna_analysis/codon_distance_2.pl \
            "${ANALYSIS_DIR}"/ath_main_transcript_stop_codon.tsv \
            >"${ANALYSIS_DIR}"/ath_main_transcript_stop_codon.yml &

    wait

    log_info "Calculating codon distances for each tissue..."
    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        perl "${NJU_SEQ_DIR}"/mrna_analysis/codon_distance_3.pl \
            "${ANALYSIS_DIR}"/ath_main_transcript_start_codon.yml \
            <(sed '1d' "${OUTPUT_DIR}"/mrna.filter."${tissue}"_RF.tsv) \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_start_codon_distance_bar.tsv \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_start_codon_distance.tsv &

        perl "${NJU_SEQ_DIR}"/mrna_analysis/codon_distance_3.pl \
            "${ANALYSIS_DIR}"/ath_main_transcript_stop_codon.yml \
            <(sed '1d' "${OUTPUT_DIR}"/mrna.filter."${tissue}"_RF.tsv) \
            "${ANALYSIS_DIR}"/mrna.filter."${tissue}"_stop_codon_distance_bar.tsv \
            >"${ANALYSIS_DIR}"/mrna.filter."${tissue}"_stop_codon_distance.tsv &
    done
    wait

    rm -f "${ANALYSIS_DIR}"/ath.gff3.tmp "${ANALYSIS_DIR}"/ath_cds.gff3.tmp

    log_info "Generating codon distance comparison plot..."
    perl "${NJU_SEQ_DIR}"/presentation/codon_distance.R \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_leaf_start_codon_distance_bar.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_root_start_codon_distance_bar.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_stem_start_codon_distance_bar.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_flower_start_codon_distance_bar.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_leaf_stop_codon_distance_bar.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_root_stop_codon_distance_bar.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_stem_stop_codon_distance_bar.tsv \
        "${ANALYSIS_DIR}"/mrna.filter.Ath_flower_stop_codon_distance_bar.tsv \
        "${FIGURES_DIR}"/codon_context_comparison.pdf

    log_info "Codon analysis completed"
}

################################################################################
# Step 9: Comparative and Statistical Analysis
################################################################################

step9_venn_analysis() {
    log_info "=== Step 9: Comparative and Statistical Analysis ==="

    perl "${NJU_SEQ_DIR}"/presentation/point_venn.sh \
        <(sed '1d' "${OUTPUT_DIR}"/mrna.filter.Ath_leaf_RF.tsv) \
        <(sed '1d' "${OUTPUT_DIR}"/mrna.filter.Ath_root_RF.tsv) \
        <(sed '1d' "${OUTPUT_DIR}"/mrna.filter.Ath_stem_RF.tsv) \
        <(sed '1d' "${OUTPUT_DIR}"/mrna.filter.Ath_flower_RF.tsv) \
        "${FIGURES_DIR}"/tissue_venn_diagram.pdf

    log_info "Venn analysis completed"
}

step9_vertical_barplot() {
    log_info "Generating vertical bar plots..."

    perl "${NJU_SEQ_DIR}"/presentation/vertical_barplot.R \
        "${STATS_DIR}"/nm_site_counts.tsv \
        "${FIGURES_DIR}"/nm_site_counts_barplot.pdf

    log_info "Bar plot generation completed"
}

step9_tissue_comparison() {
    log_info "Performing tissue-specific Nm site analysis..."

    perl "${NJU_SEQ_DIR}"/mrna_analysis/specific_point_matrix.pl \
        <(sed '1d' "${OUTPUT_DIR}"/mrna.filter.Ath_leaf_RF.tsv) \
        <(sed '1d' "${OUTPUT_DIR}"/mrna.filter.Ath_root_RF.tsv) \
        <(sed '1d' "${OUTPUT_DIR}"/mrna.filter.Ath_stem_RF.tsv) \
        <(sed '1d' "${OUTPUT_DIR}"/mrna.filter.Ath_flower_RF.tsv) \
        >"${ANALYSIS_DIR}"/tissue_specific_matrix.tsv

    perl "${NJU_SEQ_DIR}"/mrna_analysis/specific_point_withsocre.pl \
        "${ANALYSIS_DIR}"/tissue_specific_matrix.tsv \
        >"${ANALYSIS_DIR}"/tissue_specific_with_score.tsv

    log_info "Tissue comparison completed"
}

################################################################################
# Step 10: Report Generation
################################################################################

step10_generate_summary() {
    log_info "=== Step 10: Report Generation ==="

    cat > "${REPORTS_DIR}"/analysis_summary.md << 'EOF'
# Nm-site Analysis Summary Report

## Overview
This report summarizes the comprehensive Nm-site analysis results across
Arabidopsis thaliana tissues including leaf, root, stem, and flower.

## Analysis Modules

### 1. Basic Nm-site Identification
- Total Nm sites identified per tissue
- Score distribution
- Sequence context analysis

### 2. Regional Distribution
- UTR5, UTR3, CDS, and intron distribution
- Tissue-specific regional preferences

### 3. Alternative Splicing
- Splice junction analysis
- Novel splice site identification
- Differential splicing patterns

### 4. Gene Annotation
- Gene ontology enrichment
- Pathway analysis
- Transcript isoform distribution

### 5. Motif Analysis
- Sequence motif discovery
- Cross-tissue motif comparison
- Position-specific nucleotide preferences

### 6. Comparative Analysis
- Tissue-specific Nm sites
- Venn diagram comparison
- Hierarchical clustering

## Quality Metrics
- Sequencing depth
- Mapping efficiency
- Duplication rates
- Score thresholds

## Files Generated
- TSV files: Detailed numerical results
- PDF files: Visualization plots
- Statistics files: Summary statistics

EOF

    log_info "Summary report generated"
}

step10_qc_report() {
    log_info "Generating quality control report..."

    for sample in "${SAMPLES[@]}"; do
        if [[ -f "${OUTPUT_DIR}/${sample}/mrna.filter.bowtie2.log" ]]; then
            perl "${NJU_SEQ_DIR}"/tool/stat_alignment.pl \
                "${OUTPUT_DIR}"/${sample}/mrna.filter.bowtie2.log \
                >> "${STATS_DIR}"/alignment_statistics.tsv
        fi
    done

    perl "${NJU_SEQ_DIR}"/tool/draw_table.R \
        "${STATS_DIR}"/alignment_statistics.tsv \
        "${FIGURES_DIR}"/alignment_summary.pdf

    log_info "QC report generated"
}

################################################################################
# Main Execution
################################################################################

main() {
    log_info "=========================================="
    log_info "NJU-seq Extended Analysis Pipeline"
    log_info "=========================================="

    check_dependencies
    setup_directories

    log_info "Configuration:"
    log_info "  - Samples: ${#SAMPLES[@]} samples"
    log_info "  - Threads: $THREADS"
    log_info "  - Output dir: $ANALYSIS_DIR"

    step5_region_distribution
    step5_filter_by_region
    step5_exon_intron_analysis

    step6_identify_alter_splice
    step6_judge_splice_exon
    step6_splice_statistics

    step7_gene_expression
    step7_transcript_selection
    step7_filter_overlapping_genes

    step8_motif_analysis
    step8_motif_comparison
    step8_sequence_logo
    step8_codon_analysis

    step9_venn_analysis
    step9_vertical_barplot
    step9_tissue_comparison

    step10_generate_summary
    step10_qc_report

    log_info "=========================================="
    log_info "Extended analysis completed successfully!"
    log_info "Results saved to: ${ANALYSIS_DIR}"
    log_info "=========================================="
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
