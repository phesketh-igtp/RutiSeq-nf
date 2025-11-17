include { SYLPH_CLASSIFICATION }    from '../modules/control_wf/sylph/read_classification/main.nf'
include { SKA_PROFILING }           from '../modules/control_wf/ska/profiling/main.nf'
include { READS_STATS }             from '../modules/control_wf/seqkit/stats/main.nf'
include { READ_TAXONOMY_QC_REPORT } from '../modules/control_wf/report/main.nf'
include { MULTIQC }                 from '../modules/control_wf/multiqc/main.nf'

workflow CONTROL_WF {

    /*
        Define the inputs from main.nf
    */

    take:
        samplesheet
        sylph_db

    main:

        /*
            Read Classification and QC controls
        */ 

        // SYLPH classification of reads
            SYLPH_CLASSIFICATION( 
                            samplesheet,
                            sylph_db
                            )

        // SKA classification of reads
            SKA_PROFILING(samplesheet)

        // FASTQC and MULTIQC
            MULTIQC(samplesheet)

        // Read Statistics
            READS_STATS(samplesheet)

        // Report of read taxonomy and quality for controls
            READ_TAXONOMY_QC_REPORT( 
                                    SYLPH_CLASSIFICATION.out.sylph_out,
                                    SKA_PROFILING.out.ska_out,
                                    READS_STATS.out.reads_stats_out,
                                    samplesheet
                                    )

    emit:
        sylph_results   = SYLPH_CLASSIFICATION.out.sylph_out
        ska_results     = SKA_PROFILING.out.ska_out
        reads_qc_html   = READ_TAXONOMY_QC_REPORT.out.html

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