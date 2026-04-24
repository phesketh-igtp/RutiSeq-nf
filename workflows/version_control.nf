include { MTBSEQ_ENV }      from '../modules/version_control/mtbseq/main.nf'
include { READQC_EV }       from '../modules/version_control/readQC/main.nf'
include { TBPROFILER_ENV }  from '../modules/version_control/tbprofiler/main.nf'
include { PHYLOGENY_ENV }   from '../modules/version_control/phylogeny/main.nf'
include { R_STATS_ENV }     from '../modules/version_control/r_stats/main.nf'
include { R_PHYLOGENY_ENV } from '../modules/version_control/r_phylogeny/main.nf'
include { SNIPPY_ENV }      from '../modules/version_control/snippy/main.nf'

workflow VERSION_CONTROL {

    /*
        Generate the YAML files for the environments for all the containers
    */

    take:
        runID

    main:
        // Generate environment YAML files for all tools
        MTBSEQ_ENV(runID)
        READQC_EV(runID)
        TBPROFILER_ENV(runID)
        PHYLOGENY_ENV(runID)
        R_STATS_ENV(runID)
        R_PHYLOGENY_ENV(runID)
        SNIPPY_ENV(runID)

}

/*
@author: Poppy J Hesketh Best
@date: 2026-04-14
@version: 1.0.0
@description: 
    Just opens the containers and exports the yaml files with all the versions into the
        results directory
@changelog
        v1.0.0-2026-04-14: Initial version
*/