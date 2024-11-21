#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MTBSEQ_SINGLE }           from './modules/local/mtbseq/single/main'
include { TBPROFILER_DB }           from './modules/local/tbprofiler/db/main'
include { TBPROFILER_PROFILE_TBDB } from './modules/local/tbprofiler/profile.tbdb/main'
include { TBPROFILER_PROFILE_WHO }  from './modules/local/tbprofiler/profile.who/main'

workflow {
    // Create channel from sample sheet
    if (params.samplesheet == null) {
        error "Please provide a samplesheet CSV file with --samplesheet"
    }

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sampleID, file(row.forward_path), file(row.reverse_path)) }
        .set { samples_ch }

    // Run TBPROFILER_PROFILE_DB and emit a dummy value
    TBPROFILER_DB()
    db_done = TBPROFILER_DB.out.collect()

    // Run TBPROFILER_PROFILE_TBDB after TBPROFILER_PROFILE_DB is done
    TBPROFILER_PROFILE_TBDB(samples_ch.combine(db_done))

    // Prepare and run TBPROFILER_PROFILE_WHO after TBPROFILER_PROFILE_DB is done
    vcf_ch = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
    TBPROFILER_PROFILE_WHO(vcf_ch)

    // Run MTBSEQ_SINGLE (this can run independently if it doesn't depend on the DB update)
    /// MTBSEQ_SINGLE(samples_ch)

}