process COVERAGPERTARGET{

    label 'low'

    input:
    path cov_region
    path threshold_region

    output:
    path "${params.sample}_mqc.yml", emit: mqc_yml

    script:
    """
    python ${projectDir}/lib/templates/coverage.py ${cov_region} ${threshold_region}
    """
}