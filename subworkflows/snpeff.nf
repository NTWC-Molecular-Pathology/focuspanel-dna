process SNPEFF{ 

    label "low"


    input:
    path (vcf)
    path (vcf_index)

    output:
    path("*.snpEff.vcf"), emit: snpeff_vcf
    path "*.csv"                      , emit: report
    path "*.html"                     , emit: summary_html
    path "*.genes.txt"                , emit: genes_txt
    path "versions.yml"               , emit: versions
    
    script:
    def prefix = task.ext.prefix
    """
    snpEff "-Xmx${task.memory.toGiga()}g" GRCh37.p13 -d -v \
    ${vcf} \
    -nodownload \
    -csvStats ${params.sample}.csv \
    > ${prefix}.snpEff.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpEff: \$(echo \$(snpEff -version 2>&1) | sed 's/.*SnpEff //; s/ 202.*\$//')
    END_VERSIONS
    """
}