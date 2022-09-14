process SORTANDFIXTAG {

    label "low"

    input:
    path (bam)

    output:
    path "${params.sample}.bam", emit: sortandfixedbam
    path "${params.sample}.bai"
    path "${params.sample}.bam.md5"
    path  "versions.yml"           , emit: versions


    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    SortSam \
    --INPUT ${bam} \
    --OUTPUT /dev/stdout \
    --SORT_ORDER "coordinate" \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false \
    | \
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    SetNmMdAndUqTags \
    --INPUT /dev/stdin \
    --OUTPUT ${params.sample}.bam \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true \
    --REFERENCE_SEQUENCE ${params.reference}/${params.hg19}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}