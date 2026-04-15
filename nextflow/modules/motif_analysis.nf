/*
 * Motif Analysis Module
 *
 * Identifies and analyzes sequence motifs around Nm sites to
 * understand sequence preferences.
 *
 * Input:
 *   - scored: Scored Nm site files
 *   - motif_info: Motif annotation database
 *   - tissues: List of tissue types
 *   - njuseq_dir: Path to NJU_seq scripts
 *   - output_dir: Output directory
 *
 * Output:
 *   - motif: Motif analysis results
 */

process MOTIF_ANALYSIS {
    tag "MOTIF_${tissue}"
    label 'scoring'

    publishDir "${params.output_dir}/motif_analysis",
        mode: params.publish_mode,
        pattern: "*.motif.tsv",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(tissue), path(scored_file)
    path motif_info
    val tissues
    path njuseq_dir
    val output_dir

    output:
    tuple val(tissue), path("mrna.filter.${tissue}_motif.tsv"), emit: motif
    path "mrna.filter.${tissue}_motif.tsv", emit: tsv

    script:
    """
    perl ${njuseq_dir}/mrna_analysis/motif_nm.pl \
        ${motif_info} \
        < ${scored_file} \
        > mrna.filter.${tissue}_motif.tsv
    """
}
