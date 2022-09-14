process MUTECT2 {

    label "low"

    input:
    path (bam)
    path (bai)

    output:
        path ("*.vcf.gz"), emit: vcf
        path("*.vcf.gz.stats"), emit: stats
        path("*.f1r2.tar.gz"), emit: f1r2
        path("*.vcf.gz.tbi"), emit: vcf_tbi
        path("*.bam"), emit: mutect2_bam
        path("*.bai"), emit: mutect2_bai
        path  "versions.yml"           , emit: versions
    
    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: ''
    """
    gatk --java-options "-Xmx10g" \
    Mutect2 \
    -L ${params.analysis}/${params.bed_file} \
    -R ${params.reference}/${params.hg19} \
    -I ${bam} \
    -O ${prefix}.vcf.gz \
    -germline-resource ${params.reference}/${params.af_only_gnomad_vcf} \
    -bamout ${prefix}.bam \
    --f1r2-tar-gz ${prefix}.f1r2.tar.gz \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    $args \
    --callable-depth 10 \
    --create-output-variant-index true \
    --max-population-af 1.00 \
    --genotype-germline-sites true \
    --min-base-quality-score 10
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}