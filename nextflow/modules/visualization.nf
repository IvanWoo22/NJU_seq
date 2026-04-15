/*
 * Visualization Module
 *
 * Generates publication-quality plots and figures for Nm-site analysis results.
 *
 * Input:
 *   - scored: Scored Nm site files
 *   - regional: Regional distribution files
 *   - tissues: List of tissue types
 *   - output_dir: Output directory
 *
 * Output:
 *   - plots: Generated plot files
 */

process VISUALIZATION {
    tag "VIZ_${tissue}"
    label 'visualization'

    publishDir "${params.output_dir}/figures",
        mode: params.publish_mode,
        pattern: "*.pdf",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(tissue), path(scored_file)
    tuple val(tissue), path(region_file)
    val tissues
    val output_dir

    output:
    path "mrna.filter.${tissue}.signature.pdf", emit: signature_plot
    path "mrna.filter.${tissue}.region.pdf", emit: region_plot
    path "*.pdf", emit: plots

    script:
    """
    perl ${params.njuseq_dir}/presentation/signature_count.pl \
        ${scored_file} \
        mrna.filter.${tissue}.signature.pdf

    perl ${params.njuseq_dir}/presentation/point_distribution.R \
        ${region_file} \
        mrna.filter.${tissue}.region.pdf
    """
}
