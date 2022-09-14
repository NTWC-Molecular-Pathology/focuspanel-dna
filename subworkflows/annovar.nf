process ANNOVAR{
    label 'low'

    input:
    path(vcf)

    output:
    path("*.hg19_multianno.txt")  , optional:true, emit: txt
    path "versions.yml"                 , emit: versions

    script:
    def annovar = task.ext.annovar?: ''

    """
    $annovar \
    <(gunzip -c $vcf) \
    ~/annovar/humandb/ \
    -out $params.sample \
    -buildver hg19 \
    -remove \
    -protocol refGeneWithVer,gnomad211_exome,gnomad211_genome,avsnp150,exac03,cosmic96_coding,cosmic96_noncoding,dbnsfp42c,clinvar_20220320 \
    -operation g,f,f,f,f,f,f,f,f \
    -nastring . \
    -polish \
    -vcfinput

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annovar: \$($annovar --version 2>&1 |tail -n 2 | sed 's/Version: \$Date: //;s/\$//')
    END_VERSIONS
    """
}
