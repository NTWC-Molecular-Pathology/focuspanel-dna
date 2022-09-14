process BCFTOOLS_STATS {

    label 'low'

    input:
    path vcf

    output:
    path("*stats.txt"), emit: stats
    path  "versions.yml"               , emit: versions

    script:
    def prefix = task.ext.prefix ?: ''
    """
    bcftools stats $vcf > ${prefix}.bcftools_stats.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}