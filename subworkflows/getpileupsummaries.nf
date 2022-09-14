process "GETPILEUPSUMMARIES" {

    label "low" 

    input:
        path bam
        path bai


    output:
        path ("*.getpileupsummaries.table"), emit: pileup_table
        path  "versions.yml"           , emit: versions


    script:

    def prefix = task.ext.prefix ?:''

    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    GetPileupSummaries \
    -I $bam \
    -V ${params.reference}/${params.ExAc} \
    -O ${prefix}.getpileupsummaries.table \
    -L ${params.analysis}/${params.bed_file} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}