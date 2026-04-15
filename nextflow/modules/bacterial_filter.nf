/*
 * Bacterial Filter Module
 *
 * Filters bacterial alignment results to identify reads that did NOT
 * align to bacterial reference, which are considered clean reads.
 *
 * Input:
 *   - sam: Raw SAM file from bacterial alignment
 *   - njuseq_dir: Path to NJU_seq scripts
 *   - output_dir: Output directory
 *
 * Output:
 *   - filtered_list: List of clean read IDs
 */

process BACTERIAL_FILTER {
    tag "BAC_FILTER_${sample}"
    label 'filter'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "bac_rna.out.list"

    input:
    tuple val(sample), path(sam), path(log)
    path njuseq_dir
    val output_dir

    output:
    tuple val(sample), path("bac_rna.out.list"), emit: filtered_list
    path "bac_rna.out.list", emit: list

    script:
    """
    samtools view -h -F 128 ${sam} \
        | samtools view -F 16 - \
        | awk '\$6!="*"&&\$7="="&&\$4==\$8{print \$1 "\\t" \$3 "\\t" \$4 "\\t" \$6 "\\t" \$10}' \
        | perl ${njuseq_dir}/rrna_analysis/multimatch_judge.pl \
        | cut -f1 > bac_rna.out.list
    """
}
