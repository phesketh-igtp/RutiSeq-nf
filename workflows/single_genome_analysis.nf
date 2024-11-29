include { CHECK_EXISTING_OUTPUTS }    from '../modules/local/pre-wf-check/check_outputs/main.nf'  
include { MTBC_READ_QC }              from '../modules/local/pre-wf-check/mtbc-reads-qc/main.nf'
include { TBPROFILER_PROFILE_TBDB }   from '../modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }    from '../modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }             from '../modules/local/mtbseq/single/main.nf'
include { SNP_PROFILING_SINGLE }      from '../modules/local/snp-barcoding/single.profiling/main.nf'
include { SNP_FILTERING_SINGLE }      from '../modules/local/snp-barcoding/single.filtering/main.nf'
//include { SNP_BARCODING_SINGLE }      from '../modules/local/snp-barcoding/single.barcoding/main.nf'

workflow SINGLE_GENOME_ANALYSIS {

    take:
        samples_ch
        kaiju_names
        kaiju_nodes
        kaiju_fmi
        tbprofiler_db

    main:

        /*
            Opening message for workflow
        */

        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

log.info """
${color_purple}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
${color_red}Workflow: ${color_green}Single genome analysis${color_purple}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~${no_color}
"""

        /*
            Commence main workflow
        */
        
        // Report the samples part of the analysis
            log.info "${color_purple}Input samples:${no_color}"
                samples_ch.view { sampleID, forward, reverse -> 
                "${color_red}Sample: ${color_cyan}$sampleID${color_red}     Forward: ${color_cyan}$forward${color_red}  Reverse: ${color_cyan}$reverse${no_color}" 
            }

        // Run CHECK_EXISTING_OUTPUTS on all samples
            check_results_ch = CHECK_EXISTING_OUTPUTS(samples_ch)

        // Log the check results
            check_results_ch.view { sampleID, forward, reverse, all_outputs_exist ->
            "${color_red}CHECK_EXISTING_OUTPUTS results     Sample: ${color_cyan}$sampleID${color_red}     Outputs exist?: ${color_cyan}$all_outputs_exist${no_color}"
            }

        // Filter the samples based on the check results
            filtered_samples_ch = check_results_ch
                .filter { sampleID, forward, reverse, all_outputs_exist -> all_outputs_exist == 'false' }
                .map { sampleID, forward, reverse, all_outputs_exist -> tuple(sampleID, forward, reverse) }

        // Count the filtered samples and set as total_samples
            filtered_samples_ch
                .count()
                .set { total_samples }

        // View the count
            total_samples.view { count -> 
                "${color_red}Number of samples being analyzed: ${color_cyan}$count${no_color}" 
            }

        // Run MTBC_READ_QC on filtered samples
            MTBC_READ_QC(filtered_samples_ch,
                        kaiju_names,
                        kaiju_nodes,
                        kaiju_fmi
                        )

        // Explicitly capture the mtbc_reads output
            mtbc_reads_ch = MTBC_READ_QC.out.mtbc_reads

        // Collect all QC outputs into a single file
        /* This doesnt work - just captures the first output and ignores the rest, also not putting in the MTBseq perc. Will report that later
        all_qc_results = MTBC_READ_QC.out.qc_out
            .collectFile(name: 'all_samples_qc.tsv', keepHeader: true, sort: true)
            */

        // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
            TBPROFILER_PROFILE_TBDB(mtbc_reads_ch,
                                tbprofiler_db)

        // Prepare and run TBPROFILER_PROFILE_WHO
            tbdb_vcf_ch = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
            TBPROFILER_PROFILE_WHO(tbdb_vcf_ch, tbprofiler_db)
    
        // Run MTBSEQ_SINGLE
            MTBSEQ_SINGLE(mtbc_reads_ch)

        //WORK IN PROGRESS::module needs fixing!
        // Run SNP_PROFILING_SINGLE using the mpileup output
            SNP_PROFILING_SINGLE(MTBSEQ_SINGLE.out.mtbseq_mpileup)

        // Filter the SNPs based on Iñaki Comas labs methods ()
            SNP_FILTERING_SINGLE(SNP_PROFILING_SINGLE.out)

        /* WORK IN PROGRESS::module needs to be written! Barcoding BED needs generating!
        // Pre-classify genomes using SNP profiles
        snp_profiles_ch = SNP_PROFILING_SINGLE.out.snp_barcoding_individual_vcf
            .join(SNP_PROFILING_SINGLE.out.snp_barcoding_individual_vcf_index)
        SNP_BARCODING_SINGLE(snp_profiles_ch)
        // In the emit section:  
        */

        // Generate a progress log of the number of genomes that have completed the analysis
            SNP_FILTERING_SINGLE.out.mtbseq_vcf_annot
                                    .map { it -> 1 }
                                    .sum()
                                    .set { completed_samples }

        // Create progress log
            completed_samples
                    .combine(total_samples)
                    .subscribe { completed, total ->
                        log.info "${color_red}Progress: ${color_cyan}$completed ${color_red}/ ${color_cyan}$total ${color_red}samples completed${no_color}"
                    }

    emit: 
        // QC reads outputs
            //all_ qc_results                   = qc_results
        // TB-Profiler outputs
            tbprofiler_tbdb_json                = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_json
            tbprofiler_tbdb_txt                 = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_res
            tbprofiler_tbdb_vcf                 = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
            tbprofiler_who_json                 = TBPROFILER_PROFILE_WHO.out.tbprof_who_json
            tbprofiler_who_txt                  = TBPROFILER_PROFILE_WHO.out.tbprof_who_txt
        // MTBseq outputs
            mtbseq_bam                          = MTBSEQ_SINGLE.out.mtbseq_bam
            mtbseq_bam_index                    = MTBSEQ_SINGLE.out.mtbseq_bam_index
            mtbseq_bamlog                       = MTBSEQ_SINGLE.out.mtbseq_bamlog
            mtbseq_uncovered_positions          = MTBSEQ_SINGLE.out.mtbseq_uncovered_positions
            mtbseq_variant_positions            = MTBSEQ_SINGLE.out.mtbseq_variant_positions
            mtbseq_strain_classification        = MTBSEQ_SINGLE.out.mtbseq_strain_classification
            mtbseq_gatk_bam                     = MTBSEQ_SINGLE.out.mtbseq_gatk_bam
            mtbseq_gatk_bam_index               = MTBSEQ_SINGLE.out.mtbseq_gatk_bam_index
            mtbseq_gatk_bamlog                  = MTBSEQ_SINGLE.out.mtbseq_gatk_bamlog
            mtbseq_gatk_grp                     = MTBSEQ_SINGLE.out.mtbseq_gatk_grp
            mtbseq_gatk_intervals               = MTBSEQ_SINGLE.out.mtbseq_gatk_intervals
            mtbseq_mpileup                      = MTBSEQ_SINGLE.out.mtbseq_mpileup
            mtbseq_mpileuplog                   = MTBSEQ_SINGLE.out.mtbseq_mpileuplog
            mtbseq_position_table               = MTBSEQ_SINGLE.out.mtbseq_position_table
            mtbseq_mapping_variant_statistics   = MTBSEQ_SINGLE.out.mtbseq_mapping_variant_statistics
        // SNP Profiling outputs
            snp_profiling_vcf                   = SNP_PROFILING_SINGLE.out.mtbseq_vcf
            snp_profiling_vcf_index             = SNP_PROFILING_SINGLE.out.mtbseq_vcf_index
        // Uncomment the following line if you implement SNP_BARCODING_SINGLE in the future
        // snp_barcoding_results = SNP_BARCODING_SINGLE.out

}