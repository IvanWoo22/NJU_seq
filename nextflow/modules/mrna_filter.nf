/*
 * mRNA Filter Module
 *
 * Filters mRNA alignment results to identify high-quality,
 * uniquely mapping reads for downstream analysis.
 *
 * Input:
 *   - sam: mRNA alignment SAM file
 *   - njuseq_dir: Path to NJU_seq scripts
 *   - output_dir: Output directory
 *
 * Output:
 *   - filtered_sam: Filtered SAM file
 */

process MRNA_FILTER {
    tag "MRNA_FILTER_${sample}"
    label 'filter'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "mrna.filter.tmp",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(sam), path(log)
    path njuseq_dir
    val output_dir

    output:
    tuple val(sample), path("mrna.filter.tmp"), emit: filtered_sam
    path "mrna.filter.tmp", emit: tmp

    script:
    """
    samtools view -h -F 128 ${sam} \
        | samtools view -F 16 - \
        | awk '\$6!="*"&&\$7="="&&\$4==\$8{print \$1 "\\t" \$3 "\\t" \$4 "\\t" \$6 "\\t" \$10}' \
        | perl ${njuseq_dir}/mrna_analysis/matchquality_judge.pl \
        | perl ${njuseq_dir}/mrna_analysis/multimatch_judge.pl \
            > mrna.filter.tmp
    """
}
