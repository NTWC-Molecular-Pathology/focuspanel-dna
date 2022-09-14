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


// Run FASTQC
include { FASTQC                                           } from '../subworkflows/run_fastqc'
include { SAMTOOLS_STATS                                   } from '../subworkflows/samtools_stats'
include { QUALIMAP                                         } from '../subworkflows/qualimap'
include { MOSDEPTH                                         } from '../subworkflows/mosdepth'
include { COVERAGPERTARGET                                 } from '../subworkflows/coveragepertarget'

// GATK Preprocessing
include { FASTQTOUBAM                                      } from '../subworkflows/fastqtoubam'
include { SAMTOFASTQANDBWAMEM                              } from '../subworkflows/samtofastqandbwamem.nf'
include { MERGEBAMALIGNMENT                                } from '../subworkflows/mergebamalignment'
include { SORTANDFIXTAG                                    } from '../subworkflows/sortandfixtag'
include { LEFTALIGNINDELS                                  } from '../subworkflows/leftalignindels'
include { BASERECALIBRATOR                                 } from '../subworkflows/baserecalibrator'
include { APPLYBQSR                                        } from '../subworkflows/applybqsr'


// VARIANT CALLING
include { MUTECT2 as MUTECT2_NODOWNSAMPLE                  } from '../subworkflows/mutect2'
include { MUTECT2 as MUTECT2_DOWNSAMPLE                    } from '../subworkflows/mutect2'
include { FREEBAYES                                        } from '../subworkflows/freebayes'
include { TABIX_BGZIPTABIX                                 } from '../subworkflows/tabix_bgziptabix'
include { LOFREQ_INDELQUAL                                 } from '../subworkflows/lofreq_indelqual'
include { LOFREQ_CALL                                      } from '../subworkflows/lofreq_call'

// MUTECT2 FILTERING
include { LEARNREADORIENT as LEARNREADORIENT_NODOWNSAMPLE               } from '../subworkflows/learnreadorient'
include { FILTERMUTECTCALL as FILTERMUTECTCALL_NODOWNSAMPLE             } from '../subworkflows/filtermutectcall'
include { GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_NODOWNSAMPLE         } from '../subworkflows/getpileupsummaries'
include { CALCULATECONTAMINATION as CALCULATECONTAMINATION_NODOWNSAMPLE } from '../subworkflows/calculatecontamination'


include { LEARNREADORIENT as LEARNREADORIENT_DOWNSAMPLE                 } from '../subworkflows/learnreadorient'
include { FILTERMUTECTCALL as FILTERMUTECTCALL_DOWNSAMPLE               } from '../subworkflows/filtermutectcall'
include { GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_DOWNSAMPLE           } from '../subworkflows/getpileupsummaries'
include { CALCULATECONTAMINATION as CALCULATECONTAMINATION_DOWNSAMPLE   } from '../subworkflows/calculatecontamination'

// VCF Normalization
include { LEFTTRIMANDALIGN as LEFTTRIMANDALIGN_NODOWNSAMPLE} from '../subworkflows/lefttrimalign'
include { LEFTTRIMANDALIGN as LEFTTRIMANDALIGN_DOWNSAMPLE  } from '../subworkflows/lefttrimalign'
include { BCFTOOLS_NORM                                    } from '../subworkflows/bcftools_norm'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_LOFREQ            } from '../subworkflows/bcftools_norm'

// ANNOTATION - SNPEFF
include { SNPEFF as SNPEFF_NODOWNSAMPLE                    } from '../subworkflows/snpeff'
include { SNPEFF as SNPEFF_DOWNSAMPLE                      } from '../subworkflows/snpeff'
include { SNPEFF as SNPEFF_FREEBAYES                       } from '../subworkflows/snpeff'

include { SNPSIFTDBSNP as SNPSIFTDBSNP_NODOWNSAMPLE        } from '../subworkflows/snpsiftdbsnp'
include { SNPSIFTDBSNP as SNPSIFTDBSNP_DOWNSAMPLE          } from '../subworkflows/snpsiftdbsnp'
include { SNPSIFTDBSNP as SNPSIFTDBSNP_FREEBAYES           } from '../subworkflows/snpsiftdbsnp'

include { SNPEFFONEPERLINE as SNPEFFONEPERLINE_NODOWNSAMPLE} from '../subworkflows/snpeffoneperline'
include { SNPEFFONEPERLINE as SNPEFFONEPERLINE_DOWNSAMPLE  } from '../subworkflows/snpeffoneperline'
include { SNPEFFONEPERLINE as SNPEFFONEPERLINE_FREEBAYES   } from '../subworkflows/snpeffoneperline'

include { SNPEXTRACTFIELDS as SNPEXTRACTFIELDS_NODOWNSAMPLE} from '../subworkflows/snpextractfields'
include { SNPEXTRACTFIELDS as SNPEXTRACTFIELDS_DOWNSAMPLE  } from '../subworkflows/snpextractfields'
include { SNPEXTRACTFIELDS as SNPEXTRACTFIELDS_FREEBAYES   } from '../subworkflows/snpextractfields'

include { GENERATEREPORT as GENERATEREPORT_NODOWNSAMPLE    } from '../subworkflows/generatereport'
include { GENERATEREPORT as GENERATEREPORT_DOWNSAMPLE      } from '../subworkflows/generatereport'

include { COMBINEVCF                                       } from '../subworkflows/combinevcf'

include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_MUTECT2NODOWNSAMPLE           } from '../subworkflows/tabix_bgziptabix'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_MUTECT2DOWNSAMPLE             } from '../subworkflows/tabix_bgziptabix'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_FREEBAYES                     } from '../subworkflows/tabix_bgziptabix'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_LOFREQ                        } from '../subworkflows/tabix_bgziptabix'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_LOFREQ_NORM                   } from '../subworkflows/tabix_bgziptabix'

// ENSEMBL VEP
include { ENSEMBLVEP as ENSEMBLVEP_NODOWNSAMPLE             } from '../subworkflows/ensemblvep' 
include { ENSEMBLVEP as ENSEMBLVEP_DOWNSAMPLE               } from '../subworkflows/ensemblvep' 
include { ENSEMBLVEP as ENSEMBLVEP_FREEBAYES                } from '../subworkflows/ensemblvep' 
include { ENSEMBLVEP as ENSEMBLVEP_LOFREQ                   } from '../subworkflows/ensemblvep' 
include { BCFTOOLS_SPLITVEP                                 } from '../subworkflows/bcftools_split-vep'

// ANNOVAR
include {ANNOVAR as ANNOVAR_FREEBAYES                      } from '../subworkflows/annovar'

// MULTIQC
include { MULTIQC                                          } from '../subworkflows/multiqc'
include { CUSTOM_DUMPSOFTWAREVERSIONS                      } from '../subworkflows/dumpsoftwareversions'


workflow SOMATIC {

    take:
    paired_fastq

    main:

        // To gather all QC reports for MultiQC
        ch_reports  = Channel.empty()
        // To gather used softwares versions for MultiQC
        ch_versions = Channel.empty()
        // To gather all txt report from varianr caller
        ch_varianttxt = Channel.empty()


        // Step 1: QC
            // FASTQC
                FASTQC (paired_fastq)
                    ch_versions = ch_versions.mix(FASTQC.out.versions)
                    ch_reports  = ch_reports.mix(FASTQC.out.zip.collect{logs -> logs})

        
        /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        GATK Preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        // Step 1: 
            // CONVERT FASTQ TO UBAM
                FASTQTOUBAM(paired_fastq)
                    ch_versions = ch_versions.mix(FASTQTOUBAM.out.versions)

        // Step 2:
            // REVERT UBAM TO FASTQ AND MAP BY BWA-MEM
                SAMTOFASTQANDBWAMEM(FASTQTOUBAM.out.ubam)
                    ch_versions = ch_versions.mix(SAMTOFASTQANDBWAMEM.out.versions)
        
        // Step 3:
            // MERGE UBAM AND BAM
                MERGEBAMALIGNMENT(FASTQTOUBAM.out.ubam, SAMTOFASTQANDBWAMEM.out.unmergedbam)
                    ch_versions = ch_versions.mix(MERGEBAMALIGNMENT.out.versions)
        
        // Step 4:
            // SORT BAM AND FIX TAG
                SORTANDFIXTAG(MERGEBAMALIGNMENT.out.mergedbam)
                    ch_versions = ch_versions.mix(SORTANDFIXTAG.out.versions)
        
        // Step 5:
            // LEFT ALIGN INDELS
                LEFTALIGNINDELS(SORTANDFIXTAG.out.sortandfixedbam)
                    ch_versions = ch_versions.mix(LEFTALIGNINDELS.out.versions)
        
        // Step 6:
            // BASE RECALIBRATION
                BASERECALIBRATOR(LEFTALIGNINDELS.out.alignedbam)
                    ch_versions = ch_versions.mix(BASERECALIBRATOR.out.versions)
        
        // Step 7:
            // APPLY BQSR
                APPLYBQSR(LEFTALIGNINDELS.out.alignedbam, BASERECALIBRATOR.out.table)
                    ch_versions = ch_versions.mix(APPLYBQSR.out.versions)
        
        /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        End of GATK Preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        BAM QC
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        SAMTOOLS_STATS(APPLYBQSR.out.cleanbam,APPLYBQSR.out.cleanbambai)
            ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
            ch_reports  = ch_reports.mix(SAMTOOLS_STATS.out.stats.collect{stats -> stats})

        QUALIMAP(APPLYBQSR.out.cleanbam)
            ch_reports  = ch_reports.mix(QUALIMAP.out[0].collect{report -> report})
            ch_versions = ch_versions.mix(QUALIMAP.out.versions)

        MOSDEPTH(APPLYBQSR.out.cleanbam,APPLYBQSR.out.cleanbambai)
            ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
            ch_reports  = ch_reports.mix(MOSDEPTH.out.regions_txt.collect{report -> report})

        COVERAGPERTARGET(MOSDEPTH.out.regions_bed,MOSDEPTH.out.thresholds_bed)
            ch_reports  = ch_reports.mix(COVERAGPERTARGET.out.mqc_yml.collect{report -> report})




        /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Start of Variant Calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // MUTECT2 - NO DOWNSAMPLING
            MUTECT2_NODOWNSAMPLE(APPLYBQSR.out.cleanbam,APPLYBQSR.out.cleanbambai)
                ch_versions = ch_versions.mix(MUTECT2_NODOWNSAMPLE.out.versions)
        // MUTECT2 - DOWNSAMPLING
            MUTECT2_DOWNSAMPLE(APPLYBQSR.out.cleanbam,APPLYBQSR.out.cleanbambai)
                ch_versions = ch_versions.mix(MUTECT2_DOWNSAMPLE.out.versions)
        
        // GETPILEUPSUMMARIES
            GETPILEUPSUMMARIES_NODOWNSAMPLE(APPLYBQSR.out.cleanbam, APPLYBQSR.out.cleanbambai )
                ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NODOWNSAMPLE.out.versions)
            GETPILEUPSUMMARIES_DOWNSAMPLE(APPLYBQSR.out.cleanbam, APPLYBQSR.out.cleanbambai )
                ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_DOWNSAMPLE.out.versions)
        
        // CALCULATE CONTAMINATION
            CALCULATECONTAMINATION_NODOWNSAMPLE(GETPILEUPSUMMARIES_NODOWNSAMPLE.out.pileup_table)
                ch_versions = ch_versions.mix(CALCULATECONTAMINATION_NODOWNSAMPLE.out.versions)
            
            CALCULATECONTAMINATION_DOWNSAMPLE(GETPILEUPSUMMARIES_DOWNSAMPLE.out.pileup_table)
                ch_versions = ch_versions.mix(CALCULATECONTAMINATION_DOWNSAMPLE.out.versions)
        
        
        // LEARNREADORIENT
            LEARNREADORIENT_NODOWNSAMPLE(MUTECT2_NODOWNSAMPLE.out.f1r2)
                ch_versions = ch_versions.mix(LEARNREADORIENT_NODOWNSAMPLE.out.versions)

            LEARNREADORIENT_DOWNSAMPLE(MUTECT2_DOWNSAMPLE.out.f1r2)
                ch_versions = ch_versions.mix(LEARNREADORIENT_DOWNSAMPLE.out.versions)
        // FILTER MUTECT CALL
            FILTERMUTECTCALL_NODOWNSAMPLE(MUTECT2_NODOWNSAMPLE.out.vcf,MUTECT2_NODOWNSAMPLE.out.stats,MUTECT2_NODOWNSAMPLE.out.vcf_tbi,LEARNREADORIENT_NODOWNSAMPLE.out.read_orient_model,CALCULATECONTAMINATION_NODOWNSAMPLE.out.contam_table)
                ch_versions = ch_versions.mix(FILTERMUTECTCALL_NODOWNSAMPLE.out.versions)

            FILTERMUTECTCALL_DOWNSAMPLE(MUTECT2_DOWNSAMPLE.out.vcf,MUTECT2_DOWNSAMPLE.out.stats,MUTECT2_DOWNSAMPLE.out.vcf_tbi,LEARNREADORIENT_DOWNSAMPLE.out.read_orient_model,CALCULATECONTAMINATION_DOWNSAMPLE.out.contam_table)
                ch_versions = ch_versions.mix(FILTERMUTECTCALL_DOWNSAMPLE.out.versions)
        // NORMALIZE VCF
            LEFTTRIMANDALIGN_NODOWNSAMPLE(FILTERMUTECTCALL_NODOWNSAMPLE.out.vcf, FILTERMUTECTCALL_NODOWNSAMPLE.out.tbi)
                ch_versions = ch_versions.mix(LEFTTRIMANDALIGN_NODOWNSAMPLE.out.versions)

            LEFTTRIMANDALIGN_DOWNSAMPLE(FILTERMUTECTCALL_DOWNSAMPLE.out.vcf, FILTERMUTECTCALL_DOWNSAMPLE.out.tbi)
                ch_versions = ch_versions.mix(LEFTTRIMANDALIGN_DOWNSAMPLE.out.versions)
           
        // LOFREQ
            LOFREQ_INDELQUAL(APPLYBQSR.out.cleanbam,APPLYBQSR.out.cleanbambai)
                ch_versions = ch_versions.mix(LOFREQ_INDELQUAL.out.versions)
            LOFREQ_CALL(LOFREQ_INDELQUAL.out.bam)
                ch_versions = ch_versions.mix(LOFREQ_CALL.out.versions)
            TABIX_BGZIPTABIX_LOFREQ_NORM(LOFREQ_CALL.out.vcf)

            BCFTOOLS_NORM_LOFREQ(TABIX_BGZIPTABIX_LOFREQ_NORM.out.vcf_gz,TABIX_BGZIPTABIX_LOFREQ_NORM.out.vcf_tbi)
                ch_versions = ch_versions.mix(BCFTOOLS_NORM_LOFREQ.out.versions)


        // FREEBAYES
            FREEBAYES(APPLYBQSR.out.cleanbam,APPLYBQSR.out.cleanbambai)
                ch_versions = ch_versions.mix(FREEBAYES.out.versions)

            TABIX_BGZIPTABIX(FREEBAYES.out.freebayes_vcf)
        // NORMALIZE FREEBAYES VCF
            BCFTOOLS_NORM(TABIX_BGZIPTABIX.out.vcf_gz, TABIX_BGZIPTABIX.out.vcf_tbi)
                ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

        /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        End of Variant Calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */



        /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Start of Annotation
            1. SNPEFF annotate
            2. SNPSIFT add RSID
            3. Split effect into 1 per line
            4. SnpSift to extract fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        //SNPEFF - workflow
        
        // Mutect 2 No-Downsampling
            // SNPEF annotate
                SNPEFF_NODOWNSAMPLE(LEFTTRIMANDALIGN_NODOWNSAMPLE.out.norm_vcf,LEFTTRIMANDALIGN_NODOWNSAMPLE.out.norm_vcf_tbi)
                    ch_versions = ch_versions.mix(SNPEFF_NODOWNSAMPLE.out.versions)
                    ch_reports  = ch_reports.mix(SNPEFF_NODOWNSAMPLE.out.report.collect{report -> report})
            // SnpSift add rsid
                SNPSIFTDBSNP_NODOWNSAMPLE(SNPEFF_NODOWNSAMPLE.out.snpeff_vcf)
                    ch_versions = ch_versions.mix(SNPSIFTDBSNP_NODOWNSAMPLE.out.versions)
            // Split annotation    
                SNPEFFONEPERLINE_NODOWNSAMPLE(SNPSIFTDBSNP_NODOWNSAMPLE.out.vcf)
            // Extract Fields
                SNPEXTRACTFIELDS_NODOWNSAMPLE(SNPEFFONEPERLINE_NODOWNSAMPLE.out[0])
                    ch_versions = ch_versions.mix(SNPEXTRACTFIELDS_NODOWNSAMPLE.out.versions)
                    ch_varianttxt = ch_varianttxt.mix(SNPEXTRACTFIELDS_NODOWNSAMPLE.out.txt)
                //GENERATEREPORT_NODOWNSAMPLE(SNPEXTRACTFIELDS_NODOWNSAMPLE.out[0])
                //    ch_varianttxt = ch_varianttxt.mix(GENERATEREPORT_NODOWNSAMPLE.out.variant_txt_report)
            // Index the finalized VCF and publish
                TABIX_BGZIPTABIX_MUTECT2NODOWNSAMPLE(SNPSIFTDBSNP_NODOWNSAMPLE.out.vcf)


        // Mutect 2 Downsampling    
            // SNPEF annotate
                SNPEFF_DOWNSAMPLE(LEFTTRIMANDALIGN_DOWNSAMPLE.out.norm_vcf,LEFTTRIMANDALIGN_DOWNSAMPLE.out.norm_vcf_tbi)
                    ch_versions = ch_versions.mix(SNPEFF_DOWNSAMPLE.out.versions)
                    ch_reports  = ch_reports.mix(SNPEFF_DOWNSAMPLE.out.report.collect{report -> report})
            // SnpSift add rsid
                SNPSIFTDBSNP_DOWNSAMPLE(SNPEFF_DOWNSAMPLE.out.snpeff_vcf)
                    ch_versions = ch_versions.mix(SNPSIFTDBSNP_DOWNSAMPLE.out.versions)
            // Split annotation
                SNPEFFONEPERLINE_DOWNSAMPLE(SNPSIFTDBSNP_DOWNSAMPLE.out.vcf)
            // Extract Fields
                SNPEXTRACTFIELDS_DOWNSAMPLE(SNPEFFONEPERLINE_DOWNSAMPLE.out[0])
                    ch_versions = ch_versions.mix(SNPEXTRACTFIELDS_DOWNSAMPLE.out.versions)
                    ch_varianttxt = ch_varianttxt.mix(SNPEXTRACTFIELDS_DOWNSAMPLE.out.txt)
                //GENERATEREPORT_DOWNSAMPLE(SNPEXTRACTFIELDS_DOWNSAMPLE.out[0])
                    //ch_varianttxt = ch_varianttxt.mix(GENERATEREPORT_DOWNSAMPLE.out.variant_txt_report)
            
            // Index the finalized VCF and publish
                TABIX_BGZIPTABIX_MUTECT2DOWNSAMPLE(SNPSIFTDBSNP_DOWNSAMPLE.out.vcf)
        // Freebayes
            // SNPEF annotate
                SNPEFF_FREEBAYES( BCFTOOLS_NORM.out.norm_vcf, BCFTOOLS_NORM.out.norm_vcf_tbi)
                    ch_versions = ch_versions.mix(SNPEFF_FREEBAYES.out.versions)
                    ch_reports  = ch_reports.mix(SNPEFF_FREEBAYES.out.report.collect{report -> report})
            // SnpSift add rsid
                SNPSIFTDBSNP_FREEBAYES(SNPEFF_FREEBAYES.out.snpeff_vcf)
                    ch_versions = ch_versions.mix(SNPSIFTDBSNP_FREEBAYES.out.versions)
            // Split annotation
                SNPEFFONEPERLINE_FREEBAYES(SNPSIFTDBSNP_FREEBAYES.out.vcf)
            // Extract Fields
                SNPEXTRACTFIELDS_FREEBAYES(SNPEFFONEPERLINE_FREEBAYES.out[0])
                    ch_versions = ch_versions.mix(SNPEXTRACTFIELDS_FREEBAYES.out.versions)
                    ch_varianttxt = ch_varianttxt.mix(SNPEXTRACTFIELDS_FREEBAYES.out.txt)
            // Index the finalized VCF and publish
                TABIX_BGZIPTABIX_FREEBAYES(SNPSIFTDBSNP_FREEBAYES.out.vcf)
        
        
        // COMBINE the VCF from all output
            COMBINEVCF(ch_varianttxt.collect())

        // VEP

            //ENSEMBL
            ENSEMBLVEP_NODOWNSAMPLE(LEFTTRIMANDALIGN_NODOWNSAMPLE.out.norm_vcf)
                ch_versions = ch_versions.mix(ENSEMBLVEP_NODOWNSAMPLE.out.versions)

            ENSEMBLVEP_DOWNSAMPLE(LEFTTRIMANDALIGN_DOWNSAMPLE.out.norm_vcf)
                ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNSAMPLE.out.versions)

            ENSEMBLVEP_FREEBAYES(BCFTOOLS_NORM.out.norm_vcf)
                ch_versions = ch_versions.mix(ENSEMBLVEP_FREEBAYES.out.versions)

            ENSEMBLVEP_LOFREQ(BCFTOOLS_NORM_LOFREQ.out.norm_vcf)
                ch_versions = ch_versions.mix(ENSEMBLVEP_LOFREQ.out.versions)
            
            TABIX_BGZIPTABIX_LOFREQ(ENSEMBLVEP_LOFREQ.out.vcf)
            
            BCFTOOLS_SPLITVEP(TABIX_BGZIPTABIX_LOFREQ.out.vcf_gz,TABIX_BGZIPTABIX_LOFREQ.out.vcf_tbi)
                ch_versions = ch_versions.mix(BCFTOOLS_SPLITVEP.out.versions)

        // ANNOVAR
            ANNOVAR_FREEBAYES(BCFTOOLS_NORM.out.norm_vcf)
                ch_versions = ch_versions.mix(ANNOVAR_FREEBAYES.out.versions)

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        multiqc report generation
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_version_yaml = Channel.empty()
        
        CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
        ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()

        ch_multiqc_files =  Channel.empty().mix(ch_version_yaml,
                                            ch_reports.collect().ifEmpty([]))

        ch_multiqc_configs = Channel.from(ch_multiqc_config).mix(ch_multiqc_custom_config).ifEmpty([])
        MULTIQC(ch_multiqc_files.collect(), ch_multiqc_configs.collect())

        multiqc_report = MULTIQC.out.report.toList()
}