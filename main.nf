#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//include { CHECK_EXISTING_OUTPUTS }      from './modules/utilities/single/check.outputs/main.nf'
include { TBPROFILER_DB_UPDATE }        from './modules/local/tbprofiler/db/main.nf'
include { TBPROFILER_PROFILE_TBDB }     from './modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }      from './modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }               from './modules/local/mtbseq/single/main.nf'

workflow {
    // Create channel from sample sheet
    if (params.samplesheet == null) {
        error "Please provide a samplesheet CSV file with --samplesheet"
    }

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sampleID, file(row.forward_path), file(row.reverse_path)) }
        .set { samples_ch } // give the channel name

    // Run TBPROFILER_PROFILE_DB and emit a dummy value
    /// the who database needs to be downloaded seperately, so this just updates the db and generates
    /// a empty val(true) that forces the remaining tb-profiler steps to wait on hold
    TBPROFILER_DB_UPDATE()
    db_done = TBPROFILER_DB_UPDATE.out.collect()

    // Run TBPROFILER_PROFILE_TBDB after TBPROFILER_PROFILE_DB is done
    TBPROFILER_PROFILE_TBDB(samples_ch.combine(db_done))

    // Prepare and run TBPROFILER_PROFILE_WHO after TBPROFILER_PROFILE_DB is done
    vcf_ch = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
    TBPROFILER_PROFILE_WHO(vcf_ch)

    // Run MTBSEQ_SINGLE (this can run independently if it doesn't depend on the DB update)
    // Add a debug statement
    ///samples_ch.view { sample_id, forward, reverse -> 
    ///    "Debug: sample_id=${sample_id}, forward=${forward}, reverse=${reverse}" }
    MTBSEQ_SINGLE(samples_ch)

}