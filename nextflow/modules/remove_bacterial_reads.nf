/*
 * Remove Bacterial Reads Module
 *
 * Removes reads that aligned to bacterial reference from the original
 * FastQ files, creating clean datasets for downstream analysis.
 *
 * Input:
 *   - filtered_list: List of read IDs to remove
 *   - samples: List of sample IDs
 *   - data_dir: Original data directory
 *   - output_dir: Output directory
 *
 * Output:
 *   - filtered_reads: Tuple of sample, R1 and R2 filtered files
 */

process REMOVE_BACTERIAL_READS {
    tag "REMOVE_BAC_${sample}"
    label 'filter'

    publishDir "${params.data_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "*.filter.fq.gz",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(filtered_list)
    val(samples)
    val data_dir
    val output_dir

    output:
    tuple val(sample), path("R1.filter.fq.gz"), path("R2.filter.fq.gz"), emit: filtered_reads
    path "*.filter.fq.gz", emit: fastq_files

    script:
    """
    for read in R1 R2; do
        perl ${params.njuseq_dir}/tool/delete_fastq.pl \
            -n ${filtered_list} \
            -i ${data_dir}/${sample}/\${read}.origin.fq.gz \
            -o \${read}.filter.fq.gz
    done
    """
}
