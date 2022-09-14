process BCFTOOLS_SPLITVEP {
    
    label "low"

    input:
    path vcf
    path tbi

    output:
        path("*.txt"), emit: txt
        path  "versions.yml"           , emit: versions

        
    script:
    def prefix = task.ext.prefix ?: ''
    def args   = task.ext.args   ?: ''

    """
    bcftools \
    +split-vep \
    $vcf \
    -c 0-31 -d \
    | \
    bcftools query \
    -H \
    -f "${args}\\n" \
    -o ${prefix}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n 1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}