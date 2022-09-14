process FREEBAYES{

    label "high"

    input:
    path (bam)
    path (bai)

    output:
        path ("*.freebayes.vcf"), emit: freebayes_vcf
        path  "versions.yml"           , emit: versions
    
    script:

    def prefix = task.ext.prefix ?:''
    """
    freebayes \
    -t ${params.analysis}/${params.bed_file} \
    -f ${params.reference}/${params.hg19} \
    -v ${prefix}.freebayes.vcf \
    -F 0.01 \
    -C 10 \
    --pooled-continuous \
    ${bam}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/^version: v//')
    END_VERSIONS
    """

}