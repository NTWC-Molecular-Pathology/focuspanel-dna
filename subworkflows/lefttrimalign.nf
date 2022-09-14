process LEFTTRIMANDALIGN {

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
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    LeftAlignAndTrimVariants \
    --split-multi-allelics \
    -R ${params.reference}/${params.hg19} \
    -V ${vcf} \
    -O ${prefix}.norm.vcf.gz 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}