process SNPEFFONEPERLINE {

    label "low"

    input:
    path (vcf)

    output:
    path ("*.oneperline.vcf")
    
    script:
    def prefix = task.ext.prefix ?:''
    """
    cat ${vcf} \
    | \
    ${projectDir}/lib/templates/vcfeffoneperline.pl > ${prefix}.oneperline.vcf
    """
}