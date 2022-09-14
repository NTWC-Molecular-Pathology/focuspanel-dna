#!/usr/bin/env nextflow

// ENABLE DSL2
nextflow.enable.dsl=2

// PRINT VERSION AND EXIT
if(params.version){
    println """\
         ===============================================
          T M H - N G S  S O M A T I C  P I P E L I N E
         ===============================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}


// DEFINE PATHS # these are strings which are used to define input Channels,
// but they are specified here as they may be referenced in LOGGING





// PRINT STANDARD LOGGING INFO
log.info ""
log.info "         ==============================================="
log.info "          T M H - N G S  S O M A T I C  P I P E L I N E"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ===============================================" }
else {
log.info "         ===============================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         sample name  : ${params.sample}"
log.info "         base dir     : ${baseDir}"
log.info "         input dir    : ${workflow.profile.tokenize(",").contains("test") ? "-" : "${params.input}"}"
log.info "         ref dir      : ${params.reference}"
log.info "         output dir   : ${params.output}"
log.info "         mode         : ${params.SE ? "single-end" : "paired-end"}"

log.info ""
log.info "         ==============================================="
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



////////////////////
// STAGE CHANNELS //
////////////////////

/*
 *   Channels are where you define the input for the different
 *    processes which make up a pipeline. Channels indicate
 *    the flow of data, i.e. the "route" that a file will take.
 */

// STAGE BAM FILES FROM TEST PROFILE # this establishes the test data to use with -profile test
if ( workflow.profile.tokenize(",").contains("test") ){

        include { check_test_data } from './lib/functions.nf' params(readPaths: params.readPaths, singleEnd: params.SE)
        INPUT = check_test_data(params.readPaths, params.SE)

} else {

    // STAGE INPUT CHANNELS # this defines the normal input when test profile is not in use
    paired_fastq = Channel
                    .fromFilePairs("${params.input}/${params.sample}_*_{R1,R2}_*.fastq.gz", checkIfExists: true)
    snpEff_db = files("${params.analysis}/${params.snpEff_db}", checkIfExists: true)
    cov_plot_R = Channel
                    .from("${params.cov_plot_R}/coverage_plot.R")
        

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = [
                            file("$projectDir/assets/multiqc_config.yml", checkIfExists: true),
                            file("$projectDir/assets/ntwc.png", checkIfExists: true)
                            ]
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
def multiqc_report = []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SOMATIC } from './workflows/somatic'

// WORKFLOW: Run main nf-core/sarek analysis pipeline
workflow TMHNGS_SOMATIC {
    take:
    paired_fastq

    main:
    SOMATIC (paired_fastq)
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
workflow {
    TMHNGS_SOMATIC (paired_fastq)
}





//////////////////
// END PIPELINE //
//////////////////

// WORKFLOW TRACING # what to display when the pipeline finishes
// eg. with errors
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// eg. in general
workflow.onComplete {

    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${!params.debug && workflow.success ? "(cleared)" : "" }"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    // run a small clean-up script to remove "work" directory after successful completion 
    if (!params.debug && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}
