#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// local modules
include     { MTBSEQ_SINGLE }           from './modules/nf-core/mtbseq/single/main.nf'
include     { TBPROFILER_PROFILE }      from './modules/nf-core/tbprofiler/profile/main.nf'
include     { VIEW_HEAD }               from './modules/utilities/view.head.nf'

// Define fixed parameters
    params.minbqual = 5
    params.mincovf = 4
    params.mincovr = 4
    params.minphred = 20
    params.minfreq = 75
    params.unambig = 0.1
    params.window = 10

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    log.info """
    ====================================================================================================
    RutiSeq-nf:
        A reference based analysis of Illumina generated lineage Mycobacterium 
        tuberculosis genome and identified MTB lineage and genomes into transmission clusters
    ====================================================================================================
              Input samples:    ${params.sample_sheet}
    BBDD/database directory:    ${params.outdir}
                Contig file:    ${params.config}
    ====================================================================================================
    """

    // Create channel from sample sheet
    samples_ch = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(sep: ',', header: true)
        .map { row -> 
            [row.old_name, row.sampleID, file(row.forward_path), file(row.reverse_path)]
        }

    // Run processes
    TBPROFILER_PROFILE(samples_ch)
    MTBSEQ_SINGLE(samples_ch, samples_ch.map { it[1] })

    // View outputs
    VIEW_HEAD(MTBSEQ_SINGLE.out.position_variants)
    VIEW_HEAD(MTBSEQ_SINGLE.out.position_tables)

    // Publish outputs
    publish:
    TBPROFILER_PROFILE.out.csv >> 'tbprofiler'
    TBPROFILER_PROFILE.out.json >> 'tbprofiler'
    TBPROFILER_PROFILE.out.txt >> 'tbprofiler'
    TBPROFILER_PROFILE.out.vcf >> 'tbprofiler'
    MTBSEQ_SINGLE.out.called >> 'mtbseq'
    MTBSEQ_SINGLE.out.position_tables_dir >> 'mtbseq'
    MTBSEQ_SINGLE.out.classification_dir >> 'mtbseq'
    MTBSEQ_SINGLE.out.statistics_dir >> 'mtbseq'
    MTBSEQ_SINGLE.out.statistics >> 'mtbseq'
    MTBSEQ_SINGLE.out.classification >> 'mtbseq'
    MTBSEQ_SINGLE.out.position_variants >> 'mtbseq'
    MTBSEQ_SINGLE.out.position_tables >> 'mtbseq'
    }

// Output configuration
output {
    'tbprofiler' {
        path 'tbprofiler_results'
    }
    'mtbseq' {
        path 'mtbseq_results'
    }
}