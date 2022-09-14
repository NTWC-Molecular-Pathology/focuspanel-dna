process LEFTALIGNINDELS {


    label "low"

    input:
    path(bam)

    output:
    path("${params.sample}.aligned.bam"), emit: alignedbam
    path  "versions.yml"           , emit: versions

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    LeftAlignIndels \
   -R ${params.reference}/${params.hg19} \
   -I ${bam} \
   -O ${params.sample}.aligned.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}