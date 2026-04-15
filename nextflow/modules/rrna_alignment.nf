/*
 * rRNA Alignment Module
 *
 * Aligns filtered reads against rRNA reference to identify and remove
 * ribosomal RNA contamination.
 *
 * Input:
 *   - filtered_reads: Tuple of sample, R1 and R2 filtered files
 *   - index: rRNA reference index
 *   - output_dir: Output directory
 *
 * Output:
 *   - sam: SAM file with rRNA alignments
 *   - log: Alignment log file
 */

process RRNA_ALIGNMENT {
    tag "RRNA_${sample}"
    label 'rrna_align'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "rrna.filter.*",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(r1), path(r2)
    path index
    val output_dir

    output:
    tuple val(sample), path("rrna.filter.sam"), path("rrna.filter.bowtie2.log"), emit: sam
    path "rrna.filter.sam", emit: sam_file
    path "rrna.filter.bowtie2.log", emit: log

    script:
    """
    bowtie2 -p ${task.cpus} \
        -k ${params.max_multimatches} \
        -t \
        --no-unal \
        --no-mixed \
        --no-discordant \
        --end-to-end \
        -D 20 -R 3 \
        -N 0 -L 10 \
        -i S,1,0.75 \
        --np 0 \
        --xeq \
        -x ${index} \
        -1 ${r1} \
        -2 ${r2} \
        -S rrna.filter.sam \
        2>&1 | tee rrna.filter.bowtie2.log
    """
}
