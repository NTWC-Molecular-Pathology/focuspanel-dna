process GENERATEREPORT{
    
    label "low"

    input:
    path (tsv)

    output:
    path ("*.linux_pipeline_vcf.txt"), emit: variant_txt_report
    
    script:
    def prefix = task.ext.prefix ?:''
    """    
    python ${projectDir}/lib/templates/VCFtoTsv.py $tsv > ${prefix}.linux_pipeline_vcf.txt
    """

}
