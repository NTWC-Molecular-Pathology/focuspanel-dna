process ENSEMBLVEP{
    label 'low'

    input:
    path(vcf)

    output:
    path("*.vcf")  , optional:true, emit: vcf
    path("*.tsv")  , optional:true, emit: tsv
    path "versions.yml"                 , emit: versions

    script:
    def vep = task.ext.vep ?: ''
    def args = task.ext.args
    def prefix = task.ext.prefix ?:''
    def suffix = task.ext.suffix ?:''

    """
    $vep \
    $args \
    -i ${vcf} \
    --fasta ~/.vep/homo_sapiens/107_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    -o ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$($vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
