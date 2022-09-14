process VCFTOOLS {

    label 'low'

    input:
    path vcf

    output:
    path("*.TsTv.count")             , optional:true, emit: tstv_count
    path("*.TsTv.qual")              , optional:true, emit: tstv_qual
    path("*.TsTv.summary")         , optional:true, emit: tstv_summary
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    vcftools \
    --gzvcf $vcf \
    --out $prefix \
    $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}