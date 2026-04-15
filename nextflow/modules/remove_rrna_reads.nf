/*
 * Remove rRNA Reads Module
 *
 * Removes reads that aligned to rRNA from filtered datasets,
 * creating clean mRNA-enriched datasets.
 *
 * Input:
 *   - filtered_list: List of read IDs to keep
 *   - samples: List of sample IDs
 *   - data_dir: Original data directory
 *   - output_dir: Output directory
 *
 * Output:
 *   - filtered_reads: Tuple of sample, R1 and R2 filtered files
 */

process REMOVE_RRNA_READS {
    tag "REMOVE_RRNA_${sample}"
    label 'filter'

    publishDir "${params.data_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "*.filter.mrna.fq.gz",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(list)
    val samples
    val data_dir
    val output_dir

    output:
    tuple val(sample), path("R1.filter.mrna.fq.gz"), path("R2.filter.mrna.fq.gz"), emit: filtered_reads
    path "*.filter.mrna.fq.gz", emit: fastq_files

    script:
    """
    for read in R1 R2; do
        perl ${params.njuseq_dir}/tool/delete_fastq.pl \
            -n ${list} \
            -i ${data_dir}/${sample}/\${read}.filter.fq.gz \
            -o \${read}.filter.mrna.fq.gz
    done
    """
}
