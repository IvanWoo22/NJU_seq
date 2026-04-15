/*
 * mRNA Alignment Module
 *
 * Aligns rRNA-depleted reads against protein coding transcript reference
 * to identify potential Nm modification sites.
 *
 * Input:
 *   - filtered_reads: Tuple of sample, R1 and R2 filtered files
 *   - index: mRNA reference index
 *   - output_dir: Output directory
 *
 * Output:
 *   - sam: SAM file with mRNA alignments
 *   - log: Alignment log file
 */

process MRNA_ALIGNMENT {
    tag "MRNA_${sample}"
    label 'mrna_align'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "mrna.filter.*",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), path(r1), path(r2)
    path index
    val output_dir

    output:
    tuple val(sample), path("mrna.filter.sam"), path("mrna.filter.bowtie2.log"), emit: sam
    path "mrna.filter.sam", emit: sam_file
    path "mrna.filter.bowtie2.log", emit: log

    script:
    """
    bowtie2 --end-to-end \
        -p ${task.cpus} \
        -a \
        -t \
        --no-unal \
        --no-mixed \
        --no-discordant \
        -D 20 -R 3 \
        -i S,1,0.75 \
        --maxins 120 \
        -N 0 -L 10 \
        --score-min L,0.8,-0.8 \
        --rdg 8,8 --rfg 8,8 \
        --np 12 --mp 12,12 \
        --ignore-quals \
        --xeq \
        -x ${index} \
        -1 ${r1} \
        -2 ${r2} \
        -S mrna.filter.sam \
        2>&1 | tee mrna.filter.bowtie2.log
    """
}
