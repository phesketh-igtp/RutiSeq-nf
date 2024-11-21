#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MTBSEQ_SINGLE }           from './modules/local/mtbseq/single/main'
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

    // Run MTBSEQ_SINGLE
    MTBSEQ_SINGLE(samples_ch)

    // Run TBPROFILER_PROFILE_TBDB
    TBPROFILER_PROFILE_TBDB(samples_ch)

    // Prepare TBPROFILER_PROFILE_TBDB output for TBPROFILER_PROFILE_WHO
    vcf_ch = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf

    // Run TBPROFILER_PROFILE_WHO
    TBPROFILER_PROFILE_WHO(vcf_ch)

    // View outputs
    MTBSEQ_SINGLE.out.called.view()
    MTBSEQ_SINGLE.out.position_tables_dir.view()
    MTBSEQ_SINGLE.out.position_variants.view()
    MTBSEQ_SINGLE.out.position_tables.view()
    TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_bam.view()
    TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf.view()
    TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_res.view()
    TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_json.view()
    TBPROFILER_PROFILE_TBDB.out.versions.view()
}