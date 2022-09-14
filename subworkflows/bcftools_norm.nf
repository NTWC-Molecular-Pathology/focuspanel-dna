process BCFTOOLS_NORM {
    
    label "low"

    input:
    path vcf
    path tbi

    output:
        path("*.norm.vcf.gz"), emit: norm_vcf
        path("*.norm.vcf.gz.tbi"), emit: norm_vcf_tbi
        path  "versions.yml"           , emit: versions

        
    script:
    def prefix = task.ext.prefix ?:''
    """
    bcftools \
    norm \
    --threads $task.cpus \
    --atomize \
    --fasta-ref ${params.reference}/${params.hg19} \
    --multiallelics -both\
    --output ${prefix}.norm.vcf.gz \
    --output-type z \
    $vcf \
    && \
    tabix -p vcf ${prefix}.norm.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """


}