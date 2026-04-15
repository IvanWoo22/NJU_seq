/*
 * Bacterial RNA Alignment Module
 *
 * Aligns sequencing reads against bacterial RNA reference to identify
 * and remove bacterial contamination.
 *
 * Input:
 *   - sample: Sample ID
 *   - index: Bacterial reference index path
 *   - data_dir: Raw data directory
 *   - output_dir: Output directory
 *
 * Output:
 *   - sam: Raw SAM file with alignment results
 *   - log: Alignment log file
 *
 * Process:
 *   1. Bowtie2 alignment with end-to-end mode
 *   2. Multi-mapping allowed (up to 15 hits)
 *   3. Concordant and discordant pairs retained
 */

process BACTERIAL_ALIGNMENT {
    tag "BAC_${sample}"
    label 'bacterial_align'

    publishDir "${params.output_dir}/${sample}",
        mode: params.publish_mode,
        pattern: "bac_rna.*",
        saveAs: { filename -> "${filename}" }

    input:
    tuple val(sample), val(index), val(data_dir), val(output_dir)

    output:
    tuple val(sample), path("bac_rna.raw.sam"), path("bac_rna.bowtie2.log"), emit: sam
    path "bac_rna.raw.sam", emit: sam_file
    path "bac_rna.bowtie2.log", emit: log

    script:
    """
    mkdir -p ${output_dir}

    bowtie2 --end-to-end \
        -p ${task.cpus} \
        -k ${params.max_multimatches} \
        -t \
        --no-unal \
        --no-mixed \
        --no-discordant \
        -D 20 -R 3 \
        -N 0 -L 10 \
        -i S,1,0.75 \
        --maxins 120 \
        --score-min L,-0.8,-0.8 \
        --rdg 8,8 --rfg 8,8 \
        --np 12 --mp 12,12 \
        --ignore-quals \
        --xeq \
        -x ${index} \
        -1 ${data_dir}/${sample}/R1.origin.fq.gz \
        -2 ${data_dir}/${sample}/R2.origin.fq.gz \
        -S bac_rna.raw.sam \
        2>&1 | tee bac_rna.bowtie2.log
    """
}
