process TABIX_BGZIPTABIX {
    
    label "low"


    input:
    path vcf

    output:
    path '*.gz', emit: vcf_gz
    path '*.tbi', emit: vcf_tbi

    script:
    def prefix = task.ext.prefix ?:''
    """
    bgzip  --threads ${task.cpus} -c $vcf  > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz
    """

}