/*
 * mRNA Deduplication Module
 *
 * Removes duplicate reads based on mapping position and unique
 * molecular identifiers (UMIs) to reduce PCR bias.
 *
 * Input:
 *   - filtered_sam: Filtered SAM file
 *   - samples: List of sample IDs
 *   - exon_info: Exon information file
 *   - output_dir: Output directory
 *
 * Output:
 *   - dedup_sam: Deduplicated SAM file
 */

process MRNA_DEDUP {
    tag "MRNA_DEDUP_${sample}"
    label 'scoring'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "mrna.dedup.*",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(sam)
    val samples
    path exon_info
    val output_dir

    output:
    tuple val(sample), path("mrna.dedup.filter.tmp"), emit: dedup_sam
    path "mrna.dedup.filter.tmp", emit: tmp

    script:
    """
    perl ${params.njuseq_dir}/mrna_analysis/dedup.pl \
        --refstr 'Parent=transcript:' \
        --info ${exon_info} \
        < ${sam} \
        > mrna.dedup.filter.tmp
    """
}
