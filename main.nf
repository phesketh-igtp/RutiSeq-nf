#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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
        
        =========================================
        Your Workflow Name
        =========================================
        Input samples: ${params.sample_sheet}
        Output directory: ${params.outdir}
        
        =========================================
        """

    // Create channel from sample sheet/define input channel:

    def header
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(sep: ';', header: true)
        .map { row -> 
            if (header == null) {
                header = row.keySet().collect()
                println "Header: ${header}"
            }
            [row.old_name, row.sampleID, file(row.forward_path), file(row.reverse_path)]
        }
        .set { samples_ch }

    // Run TBPROFILER_PROFILE process
    TBPROFILER_PROFILE(samples_ch)

    // You can access the outputs like this:
    TBPROFILER_PROFILE.out.csv.view()
    TBPROFILER_PROFILE.out.json.view()
    TBPROFILER_PROFILE.out.txt.view()
    TBPROFILER_PROFILE.out.vcf.view()

    // Run MTBSEQ_SINGLE process
    MTBSEQ_SINGLE(
        samples_ch,
        samples_ch.map { it[1] },  // This creates the sampleID directory
    )

    // View the outputs
    MTBSEQ_SINGLE.out.called.view()
    MTBSEQ_SINGLE.out.position_tables_dir.view()
    MTBSEQ_SINGLE.out.classification_dir.view()
    MTBSEQ_SINGLE.out.statistics_dir.view()
    MTBSEQ_SINGLE.out.statistics.view()
    MTBSEQ_SINGLE.out.classification.view()
    VIEW_HEAD(MTBSEQ_SINGLE.out.position_variants)
    VIEW_HEAD(MTBSEQ_SINGLE.out.position_tables)
    MTBSEQ_SINGLE.out.versions.view()
}
