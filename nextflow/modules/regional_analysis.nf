/*
 * Regional Analysis Module
 *
 * Analyzes the genomic distribution of Nm sites across different
 * regions (UTR5, UTR3, CDS, intron).
 *
 * Input:
 *   - scored: Scored Nm site files
 *   - transcript_region: Transcript region annotation file
 *   - tissues: List of tissue types
 *   - njuseq_dir: Path to NJU_seq scripts
 *   - output_dir: Output directory
 *
 * Output:
 *   - region: Regional distribution results
 */

process REGIONAL_ANALYSIS {
    tag "REGIONAL_${tissue}"
    label 'scoring'

    publishDir "${params.output_dir}/regional_analysis",
        mode: params.publish_mode,
        pattern: "*.region.tsv",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(tissue), path(scored_file)
    path transcript_region
    val tissues
    path njuseq_dir
    val output_dir

    output:
    tuple val(tissue), path("mrna.filter.${tissue}_region.tsv"), emit: region
    path "mrna.filter.${tissue}_region.tsv", emit: tsv

    script:
    """
    perl ${njuseq_dir}/mrna_analysis/judge_region.pl \
        ${transcript_region} \
        < ${scored_file} \
        > mrna.filter.${tissue}_region.tsv
    """
}
