process LOFREQ_INDELQUAL {

    label "low"

    input:
    path (bam)
    path (bai)

    output:
        path ("*.lofreq.bam"), emit: bam
        path  "versions.yml"           , emit: versions
    
    script:
    def prefix = task.ext.prefix ?: ''
    """
    lofreq indelqual \
    $bam \
    -f ${params.reference}/${params.hg19} \
    --dindel \
    -o ${params.sample}.lofreq.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}