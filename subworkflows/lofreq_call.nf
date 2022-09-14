process LOFREQ_CALL {

    label "high"

    input:
    path (bam)

    output:
        path ("*.lofreq.vcf"), emit: vcf
        path  "versions.yml"           , emit: versions
    
    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: ''
    """
    lofreq call \
    -f ${params.reference}/${params.hg19} \
    --bed ${params.analysis}/${params.bed_file} \
    --call-indels \
    -o ${params.sample}.lofreq.vcf \
    --force-overwrite \
    --no-default-filter \
    $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}