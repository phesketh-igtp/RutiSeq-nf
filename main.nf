#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

// Module imports
include { MTBSEQ_SINGLE }          from './modules/local/mtbseq/single/main.nf'
include { TBPROFILER_PROFILE_WHO } from './modules/local/tbprofiler/profile.who/main.nf'
include { TBPROFILER_PROFILE_TBDB } from './modules/local/tbprofiler/profile.tbdb/main.nf'
include { VIEW_HEAD }              from './modules/utilities/view.head.nf'

// Define fixed parameters
params.minbqual = 5
params.mincovf = 4
params.mincovr = 4
params.minphred = 20
params.minfreq = 75
params.unambig = 0.1
params.window = 10

workflow {
    // Create channel from sample sheet
    samples_ch = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(sep: ',', header: true)
        .map { row -> 
            [row.old_name, row.sampleID, file(row.forward_path), file(row.reverse_path)]
        }

    // Run MTBSEQ_SINGLE
    MTBSEQ_SINGLE(samples_ch, samples_ch.map { it[1] })

    // Run TBPROFILER_PROFILE_TBDB
    TBPROFILER_PROFILE_TBDB(samples_ch)

    // Prepare TBPROFILER_PROFILE_TBDB output for TBPROFILER_PROFILE_WHO
    vcf_ch = TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf
        .map { vcf -> 
            def sampleID = vcf.getName().split('\\.')[0]
            return tuple(sampleID, vcf)
        }

    // Run TBPROFILER_PROFILE_WHO
    TBPROFILER_PROFILE_WHO(vcf_ch)

    // View outputs
    VIEW_HEAD(MTBSEQ_SINGLE.out.position_variants)
    VIEW_HEAD(MTBSEQ_SINGLE.out.position_tables)
    VIEW_HEAD(TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf)

    // Publish outputs
    publish:
    MTBSEQ_SINGLE.out.position_variants to: 'mtbseq_results'
    MTBSEQ_SINGLE.out.position_tables to: 'mtbseq_results'
    TBPROFILER_PROFILE_TBDB.out.tbprof_tbdb_vcf to: 'tbprofiler_results'
    TBPROFILER_PROFILE_WHO.out.tbprof_who_res to: 'tbprofiler_results'
}