//
// Subworkflow: Download databases for relevant programs
//

include { TBPROFILER_DB }   from '../modules/prepare-samples_wf/tbprofiler-db/main'
include { SYLPH_DB }        from '../modules/prepare-samples_wf/sylph-db/main'

workflow DATABASE_ACCESS_WF {
    take:
    runID // channel: value
    
    main:
    
    // Download TB-Profiler database
    TBPROFILER_DB(runID)

    // Download Sylph database  
    SYLPH_DB(runID)

    // Versions channel
    //ch_versions = Channel.empty()
    //      ch_versions = ch_versions.mix(TBPROFILER_DB.out.versions)
    //      ch_versions = ch_versions.mix(SYLPH_DB.out.versions)
    
    emit:
        tbprofiler_db = TBPROFILER_DB.out.db    // channel: path
        sylph_db      = SYLPH_DB.out.db         // channel: path
        //versions      = ch_versions             // channel: path
}