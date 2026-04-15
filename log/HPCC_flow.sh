#!/usr/bin/env bash
################################################################################
# HPCC_flow.sh - Arabidopsis thaliana RNA-seq Analysis Pipeline
#
# Description:
#   This script performs comprehensive RNA-seq analysis for Arabidopsis thaliana
#   samples, including:
#     - Bacterial RNA contamination filtering
#     - rRNA alignment and analysis
#     - mRNA alignment and deduplication
#     - Nm-site identification and scoring
#     - Result visualization
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
readonly NJU_SEQ_DIR="${PROJECT_ROOT}/NJU_seq"

readonly SAMPLES=(NJU6220 NJU6221 NJU6222 NJU6223 \
                  NJU6224 NJU6225 NJU6226 NJU6227 \
                  NJU6228 NJU6229 NJU6230 NJU6231 \
                  NJU6232 NJU6233 NJU6234 NJU6235)

readonly RNAS=(25s 18s 5-8s)
readonly READS=(R1 R2)

readonly THREADS=24
readonly BOWTIE_THREADS=22
readonly PARALLEL_JOBS=20

readonly DATA_DIR="${PROJECT_ROOT}/NJU_seq_reads"
readonly OUTPUT_DIR="${PROJECT_ROOT}/NJU_seq_output"
readonly INDEX_DIR="${PROJECT_ROOT}/index"

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
    local deps=("bowtie2" "samtools" "parallel" "perl" "pigz")
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            log_error "Required dependency '$dep' not found"
            exit 1
        fi
    done
    log_info "All dependencies satisfied"
}

check_directories() {
    log_info "Checking directory structure..."
    local dirs=("$DATA_DIR" "$OUTPUT_DIR" "$INDEX_DIR")
    for dir in "${dirs[@]}"; do
        if [[ ! -d "$dir" ]]; then
            log_error "Directory '$dir' does not exist"
            exit 1
        fi
    done
    log_info "Directory structure verified"
}

check_sample_data() {
    local sample="$1"
    for read in "${READS[@]}"; do
        local file="${DATA_DIR}/${sample}/${read}.origin.fq.gz"
        if [[ ! -f "$file" ]]; then
            log_warn "Sample ${sample} ${read} file not found: $file"
            return 1
        fi
    done
    return 0
}

################################################################################
# Step 1: Bacterial RNA Filtering
################################################################################

step1_bacterial_alignment() {
    log_info "=== Step 1: Bacterial RNA Alignment ==="
    log_info "Aligning samples against bacterial reference..."

    for sample in "${SAMPLES[@]}"; do
        log_info "Processing $sample..."

        mkdir -p "${OUTPUT_DIR}/${sample}"

        bsub -n "$THREADS" -J "${sample}_bac" "
            bowtie2 --end-to-end \
                -p ${BOWTIE_THREADS} \
                -k 15 \
                -t \
                --no-unal \
                --no-mixed \
                --no-discordant \
                -D 20 -R 3 \
                -N 0 -L 10 \
                -i S,1,0.75 \
                --maxins 120 \
                --score-min L,-0.8,-0.8 \
                --rdg 8,8 --rfg 8,8 \
                --np 12 --mp 12,12 \
                --ignore-quals \
                --xeq \
                -x ${INDEX_DIR}/bac_rna \
                -1 ${DATA_DIR}/${sample}/R1.origin.fq.gz \
                -2 ${DATA_DIR}/${sample}/R2.origin.fq.gz \
                -S ${OUTPUT_DIR}/${sample}/bac_rna.raw.sam \
                2>&1 | tee ${OUTPUT_DIR}/${sample}/bac_rna.bowtie2.log
        "
    done

    log_info "Bacterial alignment jobs submitted"
}

step1_filter_bacterial() {
    log_info "Filtering out bacterial RNA contamination..."

    bsub -n "$THREADS" -J sam_filter_bac "
        parallel -j ${PARALLEL_JOBS} --keep-order '
            samtools view -h -F 128 ${OUTPUT_DIR}/{1}/bac_rna.raw.sam \
                | samtools view -F 16 - \
                | awk '\''$6!=\"*\"&&$7==\"=\"&&$4==$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \
                | perl ${NJU_SEQ_DIR}/rrna_analysis/multimatch_judge.pl \
                | cut -f1 >${OUTPUT_DIR}/{1}/bac_rna.out.list
        ' ::: "${SAMPLES[@]}"
    "

    log_info "Bacterial filter job submitted"
}

step1_remove_bacterial_reads() {
    log_info "Removing bacterial RNA reads from data..."

    bsub -n "$THREADS" -J fastq_filter_bac "
        parallel -j ${PARALLEL_JOBS} --keep-order '
            perl ${NJU_SEQ_DIR}/tool/delete_fastq.pl \
                -n ${OUTPUT_DIR}/{1}/bac_rna.out.list \
                -i ${DATA_DIR}/{1}/{2}.origin.fq.gz \
                -o ${DATA_DIR}/{1}/{2}.filter.fq.gz
        ' ::: "${SAMPLES[@]}" ::: "${READS[@]}"
    "

    log_info "Bacterial read removal jobs submitted"
}

################################################################################
# Step 2: rRNA Analysis
################################################################################

step2_rrna_alignment() {
    log_info "=== Step 2: rRNA Alignment ==="
    log_info "Aligning samples against rRNA reference..."

    for sample in "${SAMPLES[@]}"; do
        log_info "Processing $sample for rRNA..."

        bsub -n "$THREADS" -J "${sample}_rrna" "
            bowtie2 -p ${BOWTIE_THREADS} \
                -k 15 \
                -t \
                --no-unal \
                --no-mixed \
                --no-discordant \
                --end-to-end \
                -D 20 -R 3 \
                -N 0 -L 10 \
                -i S,1,0.75 \
                --np 0 \
                --xeq \
                -x ${INDEX_DIR}/ath/ath_rrna_total \
                -1 ${DATA_DIR}/${sample}/R1.filter.fq.gz \
                -2 ${DATA_DIR}/${sample}/R2.filter.fq.gz \
                -S ${OUTPUT_DIR}/${sample}/rrna.filter.sam \
                2>&1 | tee ${OUTPUT_DIR}/${sample}/rrna.filter.bowtie2.log
        "
    done

    log_info "rRNA alignment jobs submitted"
}

step2_rrna_count() {
    log_info "Counting rRNA read ends..."

    bsub -n "$THREADS" -J rrna_count "
        parallel -j ${PARALLEL_JOBS} --keep-order '
            perl ${NJU_SEQ_DIR}/rrna_analysis/readend_count.pl \
                ${INDEX_DIR}/ath/{2}.fa \
                ${OUTPUT_DIR}/{1}/rrna.filter.tmp \
                {2} \
                >${OUTPUT_DIR}/{1}/rrna.filter.{2}.tsv
        ' ::: "${SAMPLES[@]}" ::: "${RNAS[@]}"
    "

    log_info "rRNA count job submitted"
}

step2_rrna_score() {
    log_info "Calculating rRNA scores..."

    for rna in "${RNAS[@]}"; do
        perl ${NJU_SEQ_DIR}/rrna_analysis/score.pl \
            "${OUTPUT_DIR}"/NJU62{32..35}/rrna.filter."${rna}".tsv \
            >"${OUTPUT_DIR}"/rrna.filter."${rna}".Ath_leaf_RF.tsv \
            &

        perl ${NJU_SEQ_DIR}/rrna_analysis/score.pl \
            "${OUTPUT_DIR}"/NJU62{20..23}/rrna.filter."${rna}".tsv \
            >"${OUTPUT_DIR}"/rrna.filter."${rna}".Ath_root_RF.tsv \
            &

        perl ${NJU_SEQ_DIR}/rrna_analysis/score.pl \
            "${OUTPUT_DIR}"/NJU62{24..27}/rrna.filter."${rna}".tsv \
            >"${OUTPUT_DIR}"/rrna.filter."${rna}".Ath_stem_RF.tsv \
            &

        perl ${NJU_SEQ_DIR}/rrna_analysis/score.pl \
            "${OUTPUT_DIR}"/NJU62{28..31}/rrna.filter."${rna}".tsv \
            >"${OUTPUT_DIR}"/rrna.filter."${rna}".Ath_flower_RF.tsv \
            &

        wait
    done

    log_info "rRNA scoring completed"
}

step2_filter_rrna() {
    log_info "Filtering rRNA alignment results..."

    bsub -n "$THREADS" -J sam_filter_rrna "
        parallel -j ${PARALLEL_JOBS} --keep-order '
            samtools view -h -F 128 ${OUTPUT_DIR}/{1}/rrna.filter.sam \
                | samtools view -F 16 - \
                | awk '\''$6!=\"*\"&&$7==\"=\"&&$4==$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \
                | perl ${NJU_SEQ_DIR}/rrna_analysis/matchquality_judge.pl \
                | perl ${NJU_SEQ_DIR}/rrna_analysis/multimatch_judge.pl \
                    >${OUTPUT_DIR}/{1}/rrna.filter.tmp

            cut -f1 ${OUTPUT_DIR}/{1}/rrna.filter.tmp \
                >${OUTPUT_DIR}/{1}/rrna.filter.list
        ' ::: "${SAMPLES[@]}"
    "

    log_info "rRNA filter job submitted"
}

step2_remove_rrna_reads() {
    log_info "Removing rRNA reads from data..."

    bsub -n "$THREADS" -J fastq_filter_rrna "
        parallel -j ${PARALLEL_JOBS} --keep-order '
            perl ${NJU_SEQ_DIR}/tool/delete_fastq.pl \
                -n ${OUTPUT_DIR}/{1}/rrna.filter.list \
                -i ${DATA_DIR}/{1}/{2}.filter.fq.gz \
                -o ${DATA_DIR}/{1}/{2}.filter.mrna.fq.gz
        ' ::: "${SAMPLES[@]}" ::: "${READS[@]}"
    "

    log_info "rRNA read removal jobs submitted"
}

################################################################################
# Step 3: mRNA Analysis
################################################################################

step3_mrna_alignment() {
    log_info "=== Step 3: mRNA Alignment ==="
    log_info "Aligning samples against mRNA reference..."

    for sample in "${SAMPLES[@]}"; do
        log_info "Processing $sample for mRNA..."

        bsub -n "$THREADS" -J "${sample}_mrna" "
            bowtie2 --end-to-end \
                -p ${BOWTIE_THREADS} \
                -a \
                -t \
                --no-unal \
                --no-mixed \
                --no-discordant \
                -D 20 -R 3 \
                -i S,1,0.75 \
                --maxins 120 \
                -N 0 -L 10 \
                --score-min L,0.8,-0.8 \
                --rdg 8,8 --rfg 8,8 \
                --np 12 --mp 12,12 \
                --ignore-quals \
                --xeq \
                -x ${INDEX_DIR}/ath/ath_protein_coding \
                -1 ${DATA_DIR}/${sample}/R1.filter.mrna.fq.gz \
                -2 ${DATA_DIR}/${sample}/R2.filter.mrna.fq.gz \
                -S ${OUTPUT_DIR}/${sample}/mrna.filter.sam \
                2>&1 | tee ${OUTPUT_DIR}/${sample}/mrna.filter.bowtie2.log
        "
    done

    log_info "mRNA alignment jobs submitted"
}

step3_filter_mrna() {
    log_info "Filtering mRNA alignment results..."

    bsub -n "$THREADS" -J sam_filter_mrna "
        parallel -j ${PARALLEL_JOBS} --keep-order '
            samtools view -h -F 128 ${OUTPUT_DIR}/{1}/mrna.filter.sam \
                | samtools view -F 16 - \
                | awk '\''$6!=\"*\"&&$7==\"=\"&&$4==$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \
                | perl ${NJU_SEQ_DIR}/mrna_analysis/matchquality_judge.pl \
                | perl ${NJU_SEQ_DIR}/mrna_analysis/multimatch_judge.pl \
                    >${OUTPUT_DIR}/{1}/mrna.filter.tmp
        ' ::: "${SAMPLES[@]}"
    "

    log_info "mRNA filter job submitted"
}

step3_deduplication() {
    log_info "Performing deduplication..."

    for sample in "${SAMPLES[@]}"; do
        bsub -n 1 -J "${sample}_dedup" "
            perl ${NJU_SEQ_DIR}/mrna_analysis/dedup.pl \
                --refstr 'Parent=transcript:' \
                --info ${INDEX_DIR}/ath/exon.info \
                <${OUTPUT_DIR}/${sample}/mrna.filter.tmp \
                >${OUTPUT_DIR}/${sample}/mrna.dedup.filter.tmp
        "
    done

    log_info "Deduplication jobs submitted"
}

step3_count_uniqueness() {
    log_info "Identifying unique alignment patterns..."

    bsub -n "$THREADS" -J dedup_count "
        parallel -j ${PARALLEL_JOBS} --keep-order '
            perl ${NJU_SEQ_DIR}/mrna_analysis/almostuniquematch.pl \
                ${DATA_DIR}/{1}/R1.filter.mrna.fq.gz \
                ${OUTPUT_DIR}/{1}/mrna.dedup.filter.tmp \
                ${OUTPUT_DIR}/{1}/mrna.almostuniquematch.filter.tmp

            perl ${NJU_SEQ_DIR}/mrna_analysis/count.pl \
                ${OUTPUT_DIR}/{1}/mrna.almostuniquematch.filter.tmp \
                | sort -k1,1 -k2,2n \
                >${OUTPUT_DIR}/{1}/mrna.count.filter.tmp
        ' ::: "${SAMPLES[@]}"
    "

    log_info "Count uniqueness job submitted"
}

step3_merge_samples() {
    log_info "Merging samples and adding gene information..."

    bsub -n "$THREADS" -J merge_samples "
        parallel -j ${PARALLEL_JOBS} --keep-order '
            perl ${NJU_SEQ_DIR}/mrna_analysis/merge.pl \
                --refstr 'Parent=transcript:' \
                --geneid 'AT' \
                --transid 'AT' \
                -i ${OUTPUT_DIR}/{1}/mrna.count.filter.tmp \
                --stdout \
                <${INDEX_DIR}/ath/exon.info \
                | sort -k1,1 -k2,2n \
                    >${OUTPUT_DIR}/{1}/mrna.filter.tsv
        ' ::: "${SAMPLES[@]}"
    "

    log_info "Merge job submitted"
}

################################################################################
# Step 4: Scoring and Visualization
################################################################################

step4_scoring() {
    log_info "=== Step 4: Scoring and Visualization ==="
    log_info "Calculating mRNA scores..."

    perl ${NJU_SEQ_DIR}/mrna_analysis/score.pl \
        "${OUTPUT_DIR}"/NJU622{0..3}/mrna.filter.tsv \
        >"${OUTPUT_DIR}"/mrna.filter.Ath_root_RF.tsv \
        &

    perl ${NJU_SEQ_DIR}/mrna_analysis/score.pl \
        "${OUTPUT_DIR}"/NJU622{4..7}/mrna.filter.tsv \
        >"${OUTPUT_DIR}"/mrna.filter.Ath_stem_RF.tsv \
        &

    perl ${NJU_SEQ_DIR}/mrna_analysis/score.pl \
        "${OUTPUT_DIR}"/NJU622{8..1}/mrna.filter.tsv \
        >"${OUTPUT_DIR}"/mrna.filter.Ath_flower_RF.tsv \
        &

    perl ${NJU_SEQ_DIR}/mrna_analysis/score.pl \
        "${OUTPUT_DIR}"/NJU623{2..5}/mrna.filter.tsv \
        >"${OUTPUT_DIR}"/mrna.filter.Ath_leaf_RF.tsv \
        &

    wait

    log_info "Scoring completed"
}

step4_visualization() {
    log_info "Generating signature count visualizations..."

    for tissue in Ath_leaf Ath_root Ath_stem Ath_flower; do
        perl ${NJU_SEQ_DIR}/presentation/signature_count.pl \
            "${OUTPUT_DIR}"/mrna.filter."${tissue}"_RF.tsv \
            "${OUTPUT_DIR}"/mrna.filter."${tissue}"_RF.signature.pdf
    done

    log_info "Visualization completed"
}

################################################################################
# Main Execution
################################################################################

main() {
    log_info "=========================================="
    log_info "NJU-seq Arabidopsis Analysis Pipeline"
    log_info "=========================================="

    check_dependencies
    check_directories

    log_info "Configuration:"
    log_info "  - Samples: ${#SAMPLES[@]} samples"
    log_info "  - Threads: $THREADS"
    log_info "  - Output dir: $OUTPUT_DIR"

    step1_bacterial_alignment
    step1_filter_bacterial
    step1_remove_bacterial_reads

    step2_rrna_alignment
    step2_filter_rrna
    step2_rrna_count
    step2_rrna_score
    step2_remove_rrna_reads

    step3_mrna_alignment
    step3_filter_mrna
    step3_deduplication
    step3_count_uniqueness
    step3_merge_samples

    step4_scoring
    step4_visualization

    log_info "=========================================="
    log_info "Pipeline completed successfully!"
    log_info "=========================================="
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
