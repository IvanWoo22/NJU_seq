/*
 * rRNA Scoring Module
 *
 * Calculates normalized scores for rRNA modification sites across
 * different tissue types.
 *
 * Input:
 *   - counts: List of count TSV files grouped by tissue
 *   - tissues: List of tissue types
 *   - njuseq_dir: Path to NJU_seq scripts
 *   - output_dir: Output directory
 *
 * Output:
 *   - scored: Scored rRNA modification sites
 */

process RRNA_SCORING {
    tag "RRNA_SCORING_${tissue}"
    label 'scoring'

    publishDir "${params.output_dir}",
        mode: params.publish_mode,
        pattern: "rrna.filter.*${tissue}*.tsv",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(tissue), path(counts)
    val tissues
    path njuseq_dir
    val output_dir

    output:
    tuple val(tissue), path("rrna.filter.${tissue}.tsv"), emit: scored
    path "rrna.filter.${tissue}.tsv", emit: tsv

    script:
    """
    perl ${njuseq_dir}/rrna_analysis/score.pl \
        ${counts.join(' ')} \
        > rrna.filter.${tissue}.tsv
    """
}
