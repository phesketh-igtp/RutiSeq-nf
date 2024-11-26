#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//include { CHECK_EXISTING_OUTPUTS }      from './modules/pre-wf-check/single/check.outputs/main.nf'
include { MTBC_READ_QC }                from './modules/local/pre-wf-check/mtbc-reads-qc/main.nf'
include { TBPROFILER_DB_UPDATE }        from './modules/local/tbprofiler/db/main.nf'
include { TBPROFILER_PROFILE_TBDB }     from './modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }      from './modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }               from './modules/local/mtbseq/single/main.nf'
include { CLEANUP_MTBC_READS }          from './modules/utilities/single/cleanup-mtbc-reads/main.nf'

workflow {
    // Create channel from sample sheet
    if (params.samplesheet == null) {
            error "Please provide a samplesheet CSV file with --samplesheet"
    }

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> // checks if the sample sheet is completed
            if (row.sampleID == null || row.forward_path == null || row.reverse_path == null) {
                error "Missing required column in samplesheet: ${row}"
            } // checks if the files exist
            tuple(row.sampleID, file(row.forward_path, checkIfExists: true), file(row.reverse_path, checkIfExists: true))
        }
        .set { samples_ch }

    // Run TBPROFILER_DB_UPDATE and emit a dummy value
    //TBPROFILER_DB_UPDATE()
    //db_done = TBPROFILER_DB_UPDATE.out.collect()

    // Run MTBC_READ_QC
    MTBC_READ_QC(samples_ch)

    // Collect all QC outputs into a single file
    MTBC_READ_QC.out.qc_out
        .collectFile(name: 'all_samples_qc.tsv', keepHeader: true, sort: true)
        .set { all_qc_results }
    view(all_qc_results)

    // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
    TBPROFILER_PROFILE_TBDB(MTBC_READ_QC.out.mtbc_reads)

    // Prepare and run TBPROFILER_PROFILE_WHO
    vcf_ch = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
    TBPROFILER_PROFILE_WHO(vcf_ch)

    // Run MTBSEQ_SINGLE
    MTBSEQ_SINGLE(MTBC_READ_QC.out.mtbc_reads)

    /*
    // Wait for all processes that use the MTBC reads to complete
    // Collect all outputs that use MTBC reads
    mtbc_reads_to_delete = MTBC_READ_QC.out.mtbc_reads
    
    // Create dummy channels for the other processes
    tbprofiler_done = TBPROFILER_PROFILE_TBDB.out.map { it -> it[0] }.collect()
    mtbseq_done = MTBSEQ_SINGLE.out.map { it -> it[0] }.collect()

    // Combine all conditions
    mtbc_reads_to_delete
        .combine(tbprofiler_done)
        .combine(mtbseq_done)
        .map { id, reads, _, _ -> tuple(id, reads) }
        .set { reads_for_cleanup }

    CLEANUP_MTBC_READS(reads_for_cleanup)
    */
}