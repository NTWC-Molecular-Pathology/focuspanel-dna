process FASTQTOUBAM {

    label "low"

    input:
    tuple val(sample), path(fastq)

    output:
    path("${params.sample}.unmapped.bam"), emit: ubam
    path  "versions.yml"             , emit: versions

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    FastqToSam \
    --FASTQ ${fastq[0]} \
    --FASTQ2 ${fastq[1]} \
    --OUTPUT ${params.sample}.unmapped.bam \
    --READ_GROUP_NAME ${sample} \
    --SAMPLE_NAME ${params.sample} \
    --LIBRARY_NAME ${params.rundate} \
    --PLATFORM_UNIT ${params.platform_unit} \
    --RUN_DATE ${params.rundate} \
    --PLATFORM ${params.platform_name} \
    --SEQUENCING_CENTER ${params.sequencing_center} 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}