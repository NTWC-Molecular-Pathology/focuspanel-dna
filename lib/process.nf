#!/usr/bin/env nextflow
// This file defines individual processes (separated for portability)

process FASTQC {

    publishDir "${params.output}", mode: params.publish_dir_mode, pattern: "*.{html,zip}"
    label "low" 

    input:
    tuple val(sample), path(fastq)

    output:
    path("fastqc/*.html"), emit: html
    path("fastqc/*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    script:
    """
    mkdir -p fastqc
    fastqc -o fastqc ${fastq[0]} ${fastq[1]}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}

process PairedFastQsToUnmappedBAM {

    label "low"

    input:
    tuple val(sample), path(fastq)

    output:
    path("${params.sample}.unmapped.bam"), emit: uBam
    path  "versions.yml"             , emit: versions

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    FastqToSam \
    --FASTQ ${fastq[0]} \
    --FASTQ2 ${fastq[1]} \
    --OUTPUT ${params.sample}.unmapped.bam \
    --READ_GROUP_NAME ${sample} \
    --SAMPLE_NAME ${params.sample} \
    --LIBRARY_NAME ${params.rundate} \
    --PLATFORM_UNIT ${params.platform_unit} \
    --RUN_DATE ${params.rundate} \
    --PLATFORM ${params.platform_name} \
    --SEQUENCING_CENTER ${params.sequencing_center} 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process SamToFastqAndBwaMem {

    label "low"

    input:
    path(uBam)

    output:
    path("${params.sample}.unmerged.bam"), emit: unmergedBam
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

process MergeBamAlignment {

    label "low"

    input:
    path(uBam)
    path(Bam)

    output:
    path("${params.sample}.merged.bam")
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

process SortAndFixTag {

    label "low"

    input:
    path (bam)

    output:
    path "${params.sample}.bam"
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

process LeftAlignIndels {


    label "low"

    input:
    path(bam)

    output:
    path("${params.sample}.aligned.bam")
    path  "versions.yml"           , emit: versions

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    LeftAlignIndels \
   -R ${params.reference}/${params.hg19} \
   -I ${bam} \
   -O ${params.sample}.aligned.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process BaseRecalibrator {


    label "low" 

    input:
    path (bam)

    output:
    path ("${params.sample}.recal.table"), emit: table
    path  "versions.yml"           , emit: versions

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    BaseRecalibrator \
    -I ${bam} \
    -O ${params.sample}.recal.table \
    --tmp-dir . \
    -R ${params.reference}/${params.hg19} \
    --known-sites ${params.reference}/${params.dbsnp} \
    --known-sites ${params.reference}/${params.known_site_1} \
    --known-sites ${params.reference}/${params.known_site_2} \
    --verbosity INFO
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process ApplyBQSR {

    publishDir "${params.output}/bam", mode: params.publish_dir_mode, pattern: "*.{bam,bai}"
    label "low" 
    
    input:
    path (bam)
    path (recaltable)

    output:
    path ("${params.sample}.bam")
    path ("${params.sample}.bai")
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

process SAMTOOLS_STATS {
    
    label 'low'

    input:
    path bam
    path bai


    output:
    path("*.stats"), emit: stats
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools \
        stats \
        --threads ${task.cpus-1} \
        --reference ${params.reference}/${params.hg19} \
        ${bam} \\
        > ${params.sample}.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process QUALIMAP {
    
    publishDir "${params.output}/reports/qualimap", mode: params.publish_dir_mode, pattern: "${params.sample}"
    label "low"

    input:
        path (recalbam)

    output:
        file ("*${params.sample}")
        path  "versions.yml"           , emit: versions

    script:
    """
    qualimap --java-mem-size=${task.memory.toGiga()}G \
    bamqc \
    -bam ${recalbam} \
    --paint-chromosome-limits \
    --genome-gc-distr HUMAN \
    -gff ${params.analysis}/${params.bed_file} \
    -nt ${task.cpus} \
    -outdir ${params.sample} \
    -outformat HTML
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap --version 2>&1) | sed 's/^.*QualiMap //; s/Built.*\$//')
    END_VERSIONS
    """


}

process MOSDEPTH {

    label 'low'

    input:
    path bam
    path bai


    output:
    path('*.global.dist.txt')      , emit: global_txt
    path('*.summary.txt')          , emit: summary_txt
    path('*.region.dist.txt')      , optional:true, emit: regions_txt
    path('*.per-base.d4')          , optional:true, emit: per_base_d4
    path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    path  "versions.yml"                            , emit: versions


    script:
    def args = task.ext.args ?: ''
 
    """
    mosdepth \
        --no-per-base \
        --threads $task.cpus \
        --by ${params.analysis}/${params.bed_file} \
        --threshold 500,15000 \
        $args \
        ${params.sample} \
        $bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}

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

process FREEBAYES{


    publishDir "${params.output}/variant_calling", mode: params.publish_dir_mode, pattern: "*.{vcf}"
    label "low"

    input:
    path (bam)
    path (bai)

    output:
        path ("*.freebayes.vcf"), emit: vcf
        path  "versions.yml"           , emit: versions
    
    script:

    def prefix = task.ext.prefix ?:''
    """
    freebayes \
    -t ${params.analysis}/${params.bed_file} \
    -f ${params.reference}/${params.hg19} \
    -v ${prefix}.freebayes.vcf \
    -F 0.01 \
    -C 10 \
    --pooled-continuous \
    ${bam}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/^version: v//')
    END_VERSIONS
    """

}

process SORTVCF{
    publishDir "${params.output}/variant_calling", mode: params.publish_dir_mode, pattern: "*.{vcf}"
    label "low"

    input:
        path vcf

    output:
        path ("*.sorted.vcf"), emit: vcf
        path  "versions.yml"           , emit: versions
    
    script:

    def prefix = task.ext.prefix ?:''
    """
    gatk SortVcf \
    -I $vcf \
    -O ${prefix}.sorted.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process MUTECT2 {

    publishDir "${params.output}/variant_calling", mode: params.publish_dir_mode, pattern: "*.{vcf.gz,gz.tbi}"
    label "low"

    input:
    path (bam)
    path (bai)

    output:
        path ("*.vcf.gz"), emit: vcf
        path("*.vcf.gz.stats"), emit: stats
        path("*.f1r2.tar.gz"), emit: f1r2
        path("*.vcf.gz.tbi"), emit: vcf_tbi
        path("*.bam"), emit: mutect2_bam
        path("*.bai"), emit: mutect2_bai
        path  "versions.yml"           , emit: versions
    
    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: ''
    """
    gatk --java-options "-Xmx10g" \
    Mutect2 \
    -L ${params.analysis}/${params.bed_file} \
    -R ${params.reference}/${params.hg19} \
    -I ${bam} \
    -O ${prefix}.vcf.gz \
    -bamout ${prefix}.bam \
    --f1r2-tar-gz ${prefix}.f1r2.tar.gz \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    $args \
    --callable-depth 10 \
    --create-output-variant-index true \
    --max-population-af 1.00 \
    --genotype-germline-sites true \
    --min-base-quality-score 10
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process LEARNREADORIENT {

    label "low"

    input:
    path (f1r2)
    

    output:
        path("*.read-orientation-model.tar.gz"), emit: read_orient_model
        path  "versions.yml"           , emit: versions

    script:

     def prefix = task.ext.prefix ?: ''

    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    LearnReadOrientationModel \
    -I ${f1r2} \
    -O ${prefix}.read-orientation-model.tar.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process FILTERMUTECTCALL {

    label "high"

    input:
    path (vcf)
    path (vcf_stats)
    path (vcf_index)
    path (ob_tar_gz)

    output:
        path("*.filtered.vcf.gz"), emit: vcf
        path("*.filtered.vcf.gz.tbi"), emit: tbi
        path  "versions.yml"           , emit: versions

    script:

    def prefix = task.ext.prefix ?:''
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    FilterMutectCalls \
    -R ${params.reference}/${params.hg19} \
    -V ${vcf} \
    -O ${prefix}.filtered.vcf.gz \
    -stats ${vcf_stats} \
    --ob-priors ${ob_tar_gz} \
    --min-allele-fraction 0.01 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process LEFTTRIMANDALIGN {

    label "high"

    input:
    path vcf
    path tbi

    output:
        path("*.norm.vcf.gz")
        path("*.norm.vcf.gz.tbi")
        path  "versions.yml"           , emit: versions

        
    script:
    def prefix = task.ext.prefix ?:''
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    LeftAlignAndTrimVariants \
    --split-multi-allelics \
    -R ${params.reference}/${params.hg19} \
    -V ${vcf} \
    -O ${prefix}.norm.vcf.gz 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process SNPEFF{ 

    label "low"


    input:
    path (vcf)
    path (vcf_index)

    output:
    path("*.snpEff.vcf"), emit: vcf
    path "*.csv"                      , emit: report
    path "*.html"                     , emit: summary_html
    path "*.genes.txt"                , emit: genes_txt
    path "versions.yml"               , emit: versions
    
    script:
    def prefix = task.ext.prefix
    """
    snpEff "-Xmx${task.memory.toGiga()}g" GRCh37.p13 -d -v \
    ${vcf} \
    -csvStats ${params.sample}.csv \
    > ${prefix}.snpEff.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpEff: \$(echo \$(snpEff -version 2>&1) | sed 's/.*SnpEff //; s/ 202.*\$//')
    END_VERSIONS
    """
}

process SNPSIFTDBSNP {
    publishDir "${params.output}/variant_calling", mode: params.publish_dir_mode, pattern: "*.{vcf}"
    label "low"

    input:
    path (vcf)

    output:
    path ("*.dbsnpv138.vcf")
    path  "versions.yml"           , emit: versions
    
    script:
    def prefix = task.ext.prefix ?:''
    """
    SnpSift annotate \
    ${params.reference}/${params.dbsnp} \
    ${vcf} > ${prefix}.dbsnpv138.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SnpSift: \$(echo \$(SnpSift -version 2>&1) | sed 's/^.*SnpSift version //; s/(.*\$//')
    END_VERSIONS
    """

}

process SNPEFFONEPERLINE {

    label "low"

    input:
    path (vcf)

    output:
    path ("*.oneperline.vcf.gz")
    
    script:
    def prefix = task.ext.prefix ?:''
    """
    cat ${vcf} \
    | \
    vcfeffonerperline.sh \
    | \
    bgzip -c > ${prefix}.oneperline.vcf.gz
    """
}

process ENSEMBLVEP{
     label 'low'

    input:
    path(vcf)

    output:
    path("*.ann.vcf")  , optional:true, emit: vcf
    path("*.ann.tsv")  , optional:true, emit: tsv
    path("*.ann.json") , optional:true, emit: json
    path "versions.yml"                 , emit: versions

    script:
    def vep = task.ext.vep ?: ''
    def prefix = task.ext.prefix ?:''

    """
    ${projectDir}/$vep \
    -i ${vcf} \
    --cache --refseq --check_existing --use_given_ref \
    --force_overwrite --offline \
    --fasta ~/.vep/homo_sapiens/107_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --assembly GRCh37 --tab -e \
    -o ${prefix}.ann.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(${projectDir}/$vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}

process SNPEXTRACTFIELDS {
    label "low"
    publishDir "${params.output}/variant_calling", mode: params.publish_dir_mode, pattern: "*.{txt}"

    input:
    path (vcf)

    output:
    path ("*.txt"), emit: txt
    path  "versions.yml"           , emit: versions
    
    script:
    def prefix = task.ext.prefix ?:''
    def args = task.ext.args ?:''
    """
    SnpSift extractFields -s "," -e "." ${vcf}  \
    $args \
    > ${prefix}.${workflow.manifest.name}v${workflow.manifest.version}.final_variant.txt
     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SnpSift: \$(echo \$(SnpSift -version 2>&1) | sed 's/^.*SnpSift version //; s/(.*\$//')
    END_VERSIONS
    """
}

process GENERATEREPORT{
    
    publishDir "${params.output}/variant_calling", mode: params.publish_dir_mode
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

process VCFTOOLS {

    label 'low'

    input:
    path vcf

    output:
    path("*.TsTv.count")             , optional:true, emit: tstv_count
    path("*.TsTv.qual")              , optional:true, emit: tstv_qual
    path("*.TsTv.summary")         , optional:true, emit: tstv_summary
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    """
    vcftools \
    --gzvcf $vcf \
    --out $prefix \
    $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}

process BCFTOOLS_STATS {

    label 'low'

    input:
    path vcf

    output:
    path("*stats.txt"), emit: stats
    path  "versions.yml"               , emit: versions

    script:
    def prefix = task.ext.prefix ?: ''
    """
    bcftools stats $vcf > ${prefix}.bcftools_stats.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process CUSTOM_DUMPSOFTWAREVERSIONS {
    
    label 'low'

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    script:
    template 'dumpsoftwareversions.py'
}

process MULTIQC {
    label 'low'

    input:
    path  multiqc_files, stageAs: "?/*"
    tuple path(multiqc_config), path(multiqc_logo)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    """
    multiqc \\
        --force \\
        $config \\
        $args \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

process INDEXVCF {

    publishDir "${params.output}/variant_calling", mode: params.publish_dir_mode, pattern: "*.{tbi}"
    label "low"


    input:
    path vcf

    output:
    path "*.idx", emit: tbi

    script:
    def prefix = task.ext.prefix ?:''
    """
    gatk IndexFeatureFile \
    -I $vcf
    """

}
