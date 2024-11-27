include { CHECK_EXISTING_OUTPUTS }    from '../modules/local/pre-wf-check/check_outputs/main.nf'  
include { MTBC_READ_QC }              from '../modules/local/pre-wf-check/mtbc-reads-qc/main.nf'
include { TBPROFILER_PROFILE_TBDB }   from '../modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }    from '../modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }             from '../modules/local/mtbseq/single/main.nf'
include { SNP_PROFILING_SINGLE }      from '../modules/local/snp-barcoding/single.profiling/main.nf'
include { SNP_BARCODING_SINGLE }      from '../modules/local/snp-barcoding/single.barcoding/main.nf'



workflow SINGLE_GENOME_ANALYSIS {
    
        /*
            Define input for workflow
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
        def color_reset = '\u001B[0m'

        log.info """
        ${color_purple}
        ╔════════════════════════════════════════════════════════════════════════════╗
        ║                                                                            ║
        ║  ${color_green}Sub-workflow: Single genome analysis${color_purple}                        ║
        ║                                                                            ║
        ╚════════════════════════════════════════════════════════════════════════════╝
        """

        /*
            Commence main workflow
        */

        // Check for existing outputs
            CHECK_EXISTING_OUTPUTS(samples_ch)

        // Filter samples for each module
            samples_for_tbprofiler_tbdb = CHECK_EXISTING_OUTPUTS.out
                .filter { it[3] == "0" }
                .map { it[0..2] }

            samples_for_tbprofiler_who = CHECK_EXISTING_OUTPUTS.out
                .filter { it[4] == "0" }
                .map { it[0..2] }

            samples_for_mtbseq = CHECK_EXISTING_OUTPUTS.out
                .filter { it[5] == "0" }
                .map { it[0..2] }
                
        view(CHECK_EXISTING_OUTPUTS.out.map)

        // Run MTBC_READ_QC on filtered samples
            MTBC_READ_QC(CHECK_EXISTING_OUTPUTS.out.map { it[0..2] }, 
                        kaiju_names, 
                        kaiju_nodes, 
                        kaiju_fmi)
            view(MTBC_READ_QC.out)
         /*   // Create mtbc_reads channel by joining mtbc_forward and mtbc_reverse
            mtbc_reads = MTBC_READ_QC.out.mtbc_forward
                .join(MTBC_READ_QC.out.mtbc_reverse)
                .map { sampleID, mtbc_forward, mtbc_reverse -> 
                    tuple(sampleID, mtbc_forward, mtbc_reverse) 
                }

                // View the contents of mtbc_reads - debugging
                mtbc_reads.view { sampleID, mtbc_forward, mtbc_reverse ->
                                "Sample: $sampleID\n" +
                                "Forward: $mtbc_forward (exists: ${file(mtbc_forward).exists()})\n" +
                                "Reverse: $mtbc_reverse (exists: ${file(mtbc_reverse).exists()})\n" +
                                "---"
                            }
        */
        /*    
        // Collect all QC outputs into a single file
            all_qc_results = MTBC_READ_QC.out.qc_out
                .collectFile(name: 'all_samples_qc.tsv', keepHeader: true, sort: true)

        // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
            TBPROFILER_PROFILE_TBDB(
                mtbc_reads
                    .join(samples_for_tbprofiler_tbdb, by: 0)
                    .map { it -> tuple(it[0], it[1], it[2]) }, 
                tbprofiler_db
            )

        // Prepare and run TBPROFILER_PROFILE_WHO
            vcf_ch = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
                .join(samples_for_tbprofiler_who, by: 0)
                .map { sampleID, vcf -> tuple(sampleID, vcf) }
            TBPROFILER_PROFILE_WHO(vcf_ch, tbprofiler_db)
        
        // Run MTBSEQ_SINGLE
            MTBSEQ_SINGLE(
                mtbc_reads
                    .join(samples_for_mtbseq, by: 0)
                    .map { it -> tuple(it[0], it[1], it[2]) }
            )

        /* WORK IN PROGRESS::module needs fixing!
        // Run SNP_PROFILING_SINGLE using the mpileup output
        mpileup_ch = MTBSEQ_SINGLE.out.mpileup       
        SNP_PROFILING_SINGLE(mpileup_ch)
        */

        /* WORK IN PROGRESS::module needs to be written!
        // Pre-classify genomes using SNP profiles
        snp_profiles_ch = SNP_PROFILING_SINGLE.out.snp_barcoding_individual_vcf
            .join(SNP_PROFILING_SINGLE.out.snp_barcoding_individual_vcf_index)
        SNP_BARCODING_SINGLE(snp_profiles_ch)
        // In the emit section:  
        */

    emit:
        // QC reads outputs
        qc_results = all_qc_results
        // TB-Profiler outputs
        tbprofiler_tbdb_json = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_json
        tbprofiler_tbdb_txt = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_res
        tbprofiler_tbdb_vcf = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
        tbprofiler_who_json = TBPROFILER_PROFILE_WHO.out.tbprof_who_json
        tbprofiler_who_txt = TBPROFILER_PROFILE_WHO.out.tbprof_who_txt
        // MTBseq outputs
        mtbseq_bam_dir = MTBSEQ_SINGLE.out.bam_dir
        mtbseq_bam = MTBSEQ_SINGLE.out.bam
        mtbseq_bam_index = MTBSEQ_SINGLE.out.bam_index
        mtbseq_bamlog = MTBSEQ_SINGLE.out.bamlog
        mtbseq_position_tables_dir = MTBSEQ_SINGLE.out.position_tables_dir
        mtbseq_position_tables = MTBSEQ_SINGLE.out.position_tables
        mtbseq_classification_dir = MTBSEQ_SINGLE.out.classification_dir
        mtbseq_classification = MTBSEQ_SINGLE.out.classification
        mtbseq_statistics_dir = MTBSEQ_SINGLE.out.statistics_dir
        mtbseq_statistics = MTBSEQ_SINGLE.out.statistics
        mtbseq_called_dir = MTBSEQ_SINGLE.out.called_dir
        mtbseq_position_variants = MTBSEQ_SINGLE.out.position_variants
        mtbseq_mpileup_dir = MTBSEQ_SINGLE.out.mpileup_dir
        mtbseq_mpileup = MTBSEQ_SINGLE.out.mpileup
        /*
        snp_barcoding_results = SNP_BARCODING_SINGLE.out
        */

}