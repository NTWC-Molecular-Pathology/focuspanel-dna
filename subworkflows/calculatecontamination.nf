process "CALCULATECONTAMINATION" {

    label "low" 

    input:
        path table

    output:
        path ("*.calculatecontamination.table"), emit: contam_table
        path  "versions.yml"           , emit: versions


    script:

    def prefix = task.ext.prefix ?:''

    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    CalculateContamination \
    -I $table \
    -O ${prefix}.calculatecontamination.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}