process SNPEXTRACTFIELDS {
    label "low"

    input:
    path (vcf)

    output:
    path ("*.txt"), emit: txt
    path  "versions.yml"           , emit: versions
    
    script:
    def prefix = task.ext.prefix ?:''
    def args = task.ext.args ?:''
    """
    SnpSift extractFields -s "," -e "." ${vcf}  \
    $args \
    > ${prefix}.${workflow.manifest.name}v${workflow.manifest.version}.final_variant.txt
     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SnpSift: \$(echo \$(SnpSift -version 2>&1) | sed 's/^.*SnpSift version //; s/(.*\$//')
    END_VERSIONS
    """
}