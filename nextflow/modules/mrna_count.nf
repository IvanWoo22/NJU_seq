/*
 * mRNA Count Module
 *
 * Counts uniquely mapping reads and generates read end counts
 * for Nm-site identification.
 *
 * Input:
 *   - dedup_sam: Deduplicated SAM file
 *   - samples: List of sample IDs
 *   - data_dir: Original data directory
 *   - output_dir: Output directory
 *
 * Output:
 *   - count: Read end count file
 */

process MRNA_COUNT {
    tag "MRNA_COUNT_${sample}"
    label 'scoring'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "mrna.*.tmp",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(dedup_sam)
    val samples
    val data_dir
    val output_dir

    output:
    tuple val(sample), path("mrna.count.filter.tmp"), emit: count
    path "mrna.count.filter.tmp", emit: count_file

    script:
    """
    perl ${params.njuseq_dir}/mrna_analysis/almostuniquematch.pl \
        ${data_dir}/${sample}/R1.filter.mrna.fq.gz \
        ${dedup_sam} \
        mrna.almostuniquematch.filter.tmp

    perl ${params.njuseq_dir}/mrna_analysis/count.pl \
        mrna.almostuniquematch.filter.tmp \
        | sort -k1,1 -k2,2n \
        > mrna.count.filter.tmp
    """
}
