process QUALIMAP {
    
    publishDir "${params.output}/reports/qualimap", mode: params.publish_dir_mode, pattern: "${params.sample}"
    label "low"

    input:
        path (recalbam)

    output:
        file ("*${params.sample}")
        path  "versions.yml"           , emit: versions

    script:
    """
    qualimap --java-mem-size=${task.memory.toGiga()}G \
    bamqc \
    -bam ${recalbam} \
    --paint-chromosome-limits \
    --genome-gc-distr HUMAN \
    -gff ${params.analysis}/${params.bed_file} \
    -nt ${task.cpus} \
    -outdir ${params.sample} \
    -outformat HTML
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap --version 2>&1) | sed 's/^.*QualiMap //; s/Built.*\$//')
    END_VERSIONS
    """


}