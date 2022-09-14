process FASTQC {

    publishDir "${params.output}/fastqc", mode: params.publish_dir_mode, pattern: '*.{html,zip}'
    label "low" 

    input:
    tuple val(sample), path(fastq)

    output:
    path('*.html'), emit: html
    path('*.zip') , emit: zip
    path  "versions.yml"           , emit: versions

    script:
    """
    fastqc -o . ${fastq[0]} ${fastq[1]}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}