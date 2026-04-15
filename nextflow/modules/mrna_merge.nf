/*
 * mRNA Merge Module
 *
 * Merges sample-level results and adds gene/transcript annotation
 * information.
 *
 * Input:
 *   - counts: List of count files
 *   - samples: List of sample IDs
 *   - exon_info: Exon annotation file
 *   - output_dir: Output directory
 *
 * Output:
 *   - merged: Merged TSV file with annotations
 */

process MRNA_MERGE {
    tag "MRNA_MERGE_${sample}"
    label 'scoring'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "mrna.filter.tsv",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(count)
    val samples
    path exon_info
    val output_dir

    output:
    tuple val(sample), path("mrna.filter.tsv"), emit: merged
    path "mrna.filter.tsv", emit: tsv

    script:
    """
    perl ${params.njuseq_dir}/mrna_analysis/merge.pl \
        --refstr 'Parent=transcript:' \
        --geneid 'AT' \
        --transid 'AT' \
        -i ${count} \
        --stdout \
        < ${exon_info} \
        | sort -k1,1 -k2,2n \
        > mrna.filter.tsv
    """
}
