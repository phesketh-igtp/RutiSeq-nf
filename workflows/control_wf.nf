include { SYLPH_CLASSIFICATION }    from '../modules/control_wf/sylph/read_classification/main.nf'
include { SKA_PROFILING }           from '../modules/control_wf/ska/profiling/main.nf'
include { READS_STATS }             from '../modules/control_wf/seqkit/stats/main.nf'
include { READ_TAXONOMY_QC_REPORT } from '../modules/control_wf/report/main.nf'

workflow CONTROL_WF {

    /*
        Define the inputs from main.nf
    */

    take:
        sylph_output_ch

    main:

        /*
            Opening message for workflow
        */ 

        //def purple  = '\u001B[35m'
        def green   = '\u001B[32m'
        def red     = '\u001B[31m'
        def cyan    = '\u001B[36m'
        def purple  = '\u001B[35m'
        def no_col  = '\u001B[0m'

        // SYLPH classification of reads
            SYLPH_CLASSIFICATION( params.samplesheet )

        // SKA classification of reads
            SKA_PROFILING( params.samplesheet )

        // Read Statistics
            READS_STATS( params.samplesheet )

        // Report of read taxonomy and quality for controls
            READ_TAXONOMY_QC_REPORT( 
                                    SYLPH_CLASSIFICATION.out.sylph_output_ch,
                                    SKA_PROFILING.out.ska_dist_out,
                                    READS_STATS.out.reads_stats_out
                                    )

    emit:
        sylph_results = SYLPH_CLASSIFICATION.out.sylph_output_ch
        reads_taxonomy_qc_report_out = READ_TAXONOMY_QC_REPORT.out.reads_taxonomy_qc_report_out

}

/*
    @author: Poppy J Hesketh Best
    @date: 2025-11-13
    @version: 1.0.0
    @description: 
        This is the controls inspection wf, to check the sylph classification results for controls,
            and compare them to the sample results. Ensure that the controls behave as expected.
    @changelog
        v1.0.0-2025-11-13: Initial version
*/