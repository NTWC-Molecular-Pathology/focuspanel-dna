process LEARNREADORIENT {

    label "low"

    input:
    path (f1r2)
    

    output:
        path("*.read-orientation-model.tar.gz"), emit: read_orient_model
        path  "versions.yml"           , emit: versions

    script:

     def prefix = task.ext.prefix ?: ''

    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    LearnReadOrientationModel \
    -I ${f1r2} \
    -O ${prefix}.read-orientation-model.tar.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}