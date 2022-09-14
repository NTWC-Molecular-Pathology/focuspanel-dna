process SNPSIFTDBSNP {
    
    label "low"

    input:
    path (vcf)

    output:
    path ("*.dbsnpv138.vcf"), emit: vcf
    path  "versions.yml"           , emit: versions
    
    script:
    def prefix = task.ext.prefix ?:''
    """
    SnpSift annotate \
    ${params.reference}/${params.dbsnp} \
    ${vcf} > ${prefix}.dbsnpv138.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SnpSift: \$(echo \$(SnpSift -version 2>&1) | sed 's/^.*SnpSift version //; s/(.*\$//')
    END_VERSIONS
    """

}
