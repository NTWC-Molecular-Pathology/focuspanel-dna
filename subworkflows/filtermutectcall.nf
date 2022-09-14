process FILTERMUTECTCALL {

    label "high"

    input:
    path vcf
    path vcf_stats
    path vcf_index
    path ob_tar_gz
    path contam_table

    output:
        path("*.filtered.vcf.gz"), emit: vcf
        path("*.filtered.vcf.gz.tbi"), emit: tbi
        path  "versions.yml"           , emit: versions

    script:

    def prefix = task.ext.prefix ?:''
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    FilterMutectCalls \
    -R ${params.reference}/${params.hg19} \
    -V ${vcf} \
    -O ${prefix}.filtered.vcf.gz \
    -stats ${vcf_stats} \
    --ob-priors ${ob_tar_gz} \
    --contamination-table $contam_table \
    --min-allele-fraction 0.01 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}