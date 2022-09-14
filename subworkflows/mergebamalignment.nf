process MERGEBAMALIGNMENT {

    label "low"

    input:
    path(uBam)
    path(Bam)

    output:
    path("${params.sample}.merged.bam"), emit: mergedbam
    path  "versions.yml"           , emit: versions

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
     MergeBamAlignment \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ALIGNED_BAM ${Bam} \
    --UNMAPPED_BAM ${uBam} \
    --OUTPUT ${params.sample}.merged.bam \
    --REFERENCE_SEQUENCE ${params.reference}/${params.hg19} \
    --PAIRED_RUN true \
    --SORT_ORDER "unsorted" \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --CLIP_ADAPTERS false \
    --MAX_RECORDS_IN_RAM 2000000 \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --PROGRAM_RECORD_ID "bwamem" \
    --PROGRAM_GROUP_VERSION "${params.bwa_version}" \
    --PROGRAM_GROUP_COMMAND_LINE "${params.bwa_commandline}" \
    --PROGRAM_GROUP_NAME "bwamem" \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}