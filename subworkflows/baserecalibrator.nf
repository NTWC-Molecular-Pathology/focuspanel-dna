process BASERECALIBRATOR {


    label "low" 

    input:
    path (bam)

    output:
    path ("${params.sample}.recal.table"), emit: table
    path  "versions.yml"           , emit: versions

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    BaseRecalibrator \
    -I ${bam} \
    -O ${params.sample}.recal.table \
    --tmp-dir . \
    -R ${params.reference}/${params.hg19} \
    --known-sites ${params.reference}/${params.dbsnp} \
    --known-sites ${params.reference}/${params.known_site_1} \
    --known-sites ${params.reference}/${params.known_site_2} \
    --verbosity INFO
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}