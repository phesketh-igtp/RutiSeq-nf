include { CHECK_EXISTING_OUTPUTS }    from '../modules/local/pre-wf-check/check_outputs/main.nf'  
include { MTBC_READ_QC }              from '../modules/local/pre-wf-check/mtbc-reads-qc/main.nf'
include { TBPROFILER_PROFILE_TBDB }   from '../modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }    from '../modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }             from '../modules/local/mtbseq/single/main.nf'
include { SNP_PROFILING_SINGLE }      from '../modules/local/snp-barcoding/single.profiling/main.nf'
include { SNP_BARCODING_SINGLE }      from '../modules/local/snp-barcoding/single.barcoding/main.nf'

workflow SINGLE_GENOME_ANALYSIS {
    take:
        samples_ch
        kaiju_names
        kaiju_nodes
        kaiju_fmi
        tbprofiler_db

    main:
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

        // Run MTBC_READ_QC on filtered samples
        MTBC_READ_QC(CHECK_EXISTING_OUTPUTS.out.map { it[0..2] }, 
                    kaiju_names, 
                    kaiju_nodes, 
                    kaiju_fmi)

        // Collect all QC outputs into a single file
        all_qc_results = MTBC_READ_QC.out.qc_out
            .collectFile(name: 'all_samples_qc.tsv', keepHeader: true, sort: true)

        // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
        TBPROFILER_PROFILE_TBDB(MTBC_READ_QC.out.mtbc_reads
                                .join(samples_for_tbprofiler_tbdb, by: 0), tbprofiler_db)

        // Prepare and run TBPROFILER_PROFILE_WHO
        vcf_ch = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
            .join(samples_for_tbprofiler_who, by: 0)  // Join with samples that need WHO profiling
            .map { sampleID, vcf, forward, reverse -> tuple(sampleID, vcf) }  // Keep only sampleID and VCF
        TBPROFILER_PROFILE_WHO(vcf_ch, tbprofiler_db)
        
        // Run MTBSEQ_SINGLE
        MTBSEQ_SINGLE(MTBC_READ_QC.out.mtbc_reads
            .join(samples_for_mtbseq, by: 0))

        // Run SNP_PROFILING_SINGLE using the mpileup output
        mpileup_ch = MTBSEQ_SINGLE.out.mpileup       
        SNP_PROFILING_SINGLE(mpileup_ch)

        // Pre-classify genomes using SNP profiles
        //snp_profiles_ch = SNP_PROFILING_SINGLE.out.snp_barcoding_individual_vcf
        //    .join(SNP_PROFILING_SINGLE.out.snp_barcoding_individual_vcf_index)
        //SNP_BARCODING_SINGLE(snp_profiles_ch)

    emit:
        qc_results = all_qc_results
        tbprofiler_tbdb_json = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_json
        tbprofiler_tbdb_txt = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_res
        tbprofiler_tbdb_vcf = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
        tbprofiler_who_results = TBPROFILER_PROFILE_WHO.out
        mtbseq_results = MTBSEQ_SINGLE.out
        snp_profiling_results = SNP_PROFILING_SINGLE.out

}