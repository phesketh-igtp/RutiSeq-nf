include { TBPROFILER_COMPILE_TBDB }   from '../modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_COMPILE_WHO }    from '../modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_LINEAGE_SPLITTING }  from '../modules/local/mtbseq/single/main.nf'
include { MTBSEQ_LINEAGE_PAIRWISE }   from '../modules/local/mtbseq/single/main.nf'


workflow PAIRIWISE_GENOME_ANALYSIS {
    
    take:
        // TB-Profiler files
            tbprofiler_tbdb_json
            tbprofiler_tbdb_txt
            tbprofiler_tbdb_vcf
            tbprofiler_who_json
            tbprofiler_who_txt
        // MTBSeq files
            mtbseq_variant_positions
            mtbseq_strain_classification
            mtbseq_position_table
            mtbseq_mapping_variant_statistics

    main:

        /*
            Opening message for workflow
        */

        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

        log.info """
        ${color_purple}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ${color_red}Workflow: ${color_green}Comparative/Pairwise genome analysis${color_purple}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~${no_color}
        """

        // Compile TB-Profiler results
        TBPROFILER_COMPILE_TBDB(
                                tbprofiler_tbdb_txt,
                                tbprofiler_tbdb_json
                                )
        
        TBPROFILER_COMPILE_WHO(
                                tbprofiler_who_json
                                tbprofiler_who_txt
                                )

        MTBSEQ_LINEAGE_SPLITTING(
                                mtbseq_strain_classification,
                                mtbseq_strain_classification
                                )

    emit:
        

}