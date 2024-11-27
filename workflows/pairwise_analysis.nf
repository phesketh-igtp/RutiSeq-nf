include { CHECK_EXISTING_OUTPUTS }    from '../modules/local/pre-wf-check/check_outputs/main.nf'  
include { MTBC_READ_QC }              from '../modules/local/pre-wf-check/mtbc-reads-qc/main.nf'
include { TBPROFILER_PROFILE_TBDB }   from '../modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }    from '../modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }             from '../modules/local/mtbseq/single/main.nf'
include { SNP_PROFILING_SINGLE }      from '../modules/local/snp-barcoding/single.profiling/main.nf'
include { SNP_BARCODING_SINGLE }      from '../modules/local/snp-barcoding/single.barcoding/main.nf'

workflow PAIRIWISE_GENOME_ANALYSIS {
    
            /*
            Opening message for workflow
        */

        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_reset = '\u001B[0m'

        log.info """
        ${color_purple}
        ╔════════════════════════════════════════════════════════════════════════════╗
        ║                                                                            ║
        ║  ${color_green}Sub-workflow: Pairwise genome analysis${color_purple}                      ║
        ║                                                                            ║
        ╚════════════════════════════════════════════════════════════════════════════╝
        """

    take:

    main:

    emit:


}