process COMBINEVCF{

    publishDir "${params.output}/variant_calling", mode: params.publish_dir_mode
    label "low"

    input:
    path vcf_txt_files, stageAs: "*"

    output:
    path ("*.combine_vcf.xlsx")
    
    script:

    """    
    python ${projectDir}/lib/templates/combine_vcf.py ${params.sample}
    """

}
