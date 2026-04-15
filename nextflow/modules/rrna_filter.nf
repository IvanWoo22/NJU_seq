/*
 * rRNA Filter Module
 *
 * Filters rRNA alignment results to identify uniquely mapping reads
 * and generate a list of non-rRNA reads.
 *
 * Input:
 *   - sam: rRNA alignment SAM file
 *   - njuseq_dir: Path to NJU_seq scripts
 *   - output_dir: Output directory
 *
 * Output:
 *   - filtered_list: List of non-rRNA read IDs
 */

process RRNA_FILTER {
    tag "RRNA_FILTER_${sample}"
    label 'filter'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "rrna.filter.*",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(sam), path(log)
    path njuseq_dir
    val output_dir

    output:
    tuple val(sample), path("rrna.filter.list"), path("rrna.filter.tmp"), emit: filtered_list
    path "rrna.filter.list", emit: list
    path "rrna.filter.tmp", emit: tmp

    script:
    """
    samtools view -h -F 128 ${sam} \
        | samtools view -F 16 - \
        | awk '\$6!="*"&&\$7="="&&\$4==\$8{print \$1 "\\t" \$3 "\\t" \$4 "\\t" \$6 "\\t" \$10}' \
        | perl ${njuseq_dir}/rrna_analysis/matchquality_judge.pl \
        | perl ${njuseq_dir}/rrna_analysis/multimatch_judge.pl \
            > rrna.filter.tmp

    cut -f1 rrna.filter.tmp > rrna.filter.list
    """
}
