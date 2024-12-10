//include { CHECK_EXISTING_OUTPUTS }    from '../modules/local/pre-wf-check/check_outputs/main.nf'  
include { MTBC_READ_QC }              from '../modules/local/pre-wf-check/mtbc-reads-qc/main.nf'
include { COMBINE_QC_RESULTS }        from '../modules/local/pre-wf-check/combine-qc-results/main.nf'
include { TBPROFILER_PROFILE_TBDB }   from '../modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }    from '../modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }             from '../modules/local/mtbseq/single/main.nf'
include { SNP_PROFILING_SINGLE }      from '../modules/local/snp-barcoding/single.profiling/main.nf'
include { SNP_ANNOTATING_SINGLE }     from '../modules/local/snp-barcoding/single.annotating/main.nf'
//include { SNP_BARCODING_SINGLE }      from '../modules/local/snp-barcoding/single.barcoding/main.nf'

workflow SINGLE_WORKFLOW {

        /*
            Define the inputs from main.nf
        */

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
        2024-12-09
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~${no_color}
        """          

        // Run MTBC_READ_QC on filtered samples
            MTBC_READ_QC(samples_ch,
                        kaiju_names,
                        kaiju_nodes,
                        kaiju_fmi
                        )

        // Collect all QC results
            all_qc_results = MTBC_READ_QC.out.qc_results.map { it[1] }.collect()

        // Combine QC results : NEEDS REPAIRING
            COMBINE_QC_RESULTS( // produce TSV of read QC results
                            all_qc_results, 
                            params.runID)

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

        /* WORK IN PROGRESS::module needs to be written! Barcoding BED needs generating!
        
        // Filter the SNPs based on Iñaki Comas labs methods ()
            //SNP_ANNOTATING_SINGLE(SNP_PROFILING_SINGLE.out)

        // Pre-classify genomes using SNP profiles
            snp_profiles_ch = SNP_PROFILING_SINGLE.out.snp_barcoding_individual_vcf
                .join(SNP_PROFILING_SINGLE.out.snp_barcoding_individual_vcf_index)
            SNP_BARCODING_SINGLE(snp_profiles_ch)

        // Generate summary of paths
            WF_HANDOVER(
                        SNP_PROFILING_SINGLE.out,
                        MTBSEQ_SINGLE.out
                        TBPROFILER_PROFILE_WHO.out,
                        TBPROFILER_PROFILE_TBDB.out)

        */

    emit: 
        // QC reads outputs
            //all_ qc_results                   = qc_results
        // TB-Profiler outputs
        workflow_outputs = tuple(
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
            )
            
}