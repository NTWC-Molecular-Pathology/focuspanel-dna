process APPLYBQSR {

    publishDir "${params.output}/bam", mode: params.publish_dir_mode, pattern: "*.{bam,bai}"
    label "low" 
    
    input:
    path (bam)
    path (recaltable)

    output:
    path ("${params.sample}.bam"), emit: cleanbam
    path ("${params.sample}.bai"), emit: cleanbambai
    path  "versions.yml"           , emit: versions

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    ApplyBQSR \
    -R ${params.reference}/${params.hg19} \
    --input ${bam} \
    --output ${params.sample}.bam \
    --bqsr-recal-file ${recaltable} \
    --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --create-output-bam-md5 \
    --use-original-qualities
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
