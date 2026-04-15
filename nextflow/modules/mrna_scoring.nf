/*
 * mRNA Scoring Module
 *
 * Calculates statistical scores for Nm-site identification by comparing
 * treatment and control samples.
 *
 * Input:
 *   - merged: List of merged TSV files
 *   - tissues: List of tissue types
 *   - njuseq_dir: Path to NJU_seq scripts
 *   - output_dir: Output directory
 *
 * Output:
 *   - scored: Scored Nm sites
 */

process MRNA_SCORING {
    tag "MRNA_SCORING_${tissue}"
    label 'scoring'

    publishDir "${params.output_dir}",
        mode: params.publish_mode,
        pattern: "mrna.filter.*${tissue}*.tsv",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(tissue), path(merged_files)
    val tissues
    path njuseq_dir
    val output_dir

    output:
    tuple val(tissue), path("mrna.filter.${tissue}.tsv"), emit: scored
    path "mrna.filter.${tissue}.tsv", emit: tsv

    script:
    """
    perl ${njuseq_dir}/mrna_analysis/score.pl \
        ${merged_files.join(' ')} \
        > mrna.filter.${tissue}.tsv
    """
}
