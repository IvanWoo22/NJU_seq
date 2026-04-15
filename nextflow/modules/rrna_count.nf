/*
 * rRNA Count Module
 *
 * Counts read ends mapping to different rRNA types (25s, 18s, 5-8s).
 *
 * Input:
 *   - filtered_list: List of non-rRNA read IDs
 *   - samples: List of sample IDs
 *   - rrna_fasta: Path to rRNA FASTA files
 *   - output_dir: Output directory
 *
 * Output:
 *   - count: TSV file with read end counts
 */

process RRNA_COUNT {
    tag "RRNA_COUNT_${sample}_${rna_type}"
    label 'scoring'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "rrna.filter.*.tsv",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(list), path(tmp)
    val samples
    path rrna_fasta
    val output_dir

    output:
    tuple val(sample), val(rna_type), path("rrna.filter.${rna_type}.tsv"), emit: count
    path "rrna.filter.${rna_type}.tsv", emit: tsv

    script:
    """
    perl ${params.njuseq_dir}/rrna_analysis/readend_count.pl \
        ${rrna_fasta} \
        ${tmp} \
        ${rna_type} \
        > rrna.filter.${rna_type}.tsv
    """
}
