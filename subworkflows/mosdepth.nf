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