process SAMTOFASTQANDBWAMEM {

    label "low"

    input:
    path(uBam)

    output:
    path("${params.sample}.unmerged.bam"), emit: unmergedbam
    path("${params.sample}.bwa.stderr.log")
    path  "versions.yml"           , emit: versions

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    SamToFastq \
    --INPUT ${uBam} \
    --FASTQ /dev/stdout \
    --INTERLEAVE true \
    --NON_PF true \
    | \
    bwa mem -K 100000000 -p -v 3 -t 16 -Y ${params.reference}/${params.hg19} \
    /dev/stdin -  2> >(tee ${params.sample}.bwa.stderr.log >&2) \
    | \
    samtools view -1 - > ${params.sample}.unmerged.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}