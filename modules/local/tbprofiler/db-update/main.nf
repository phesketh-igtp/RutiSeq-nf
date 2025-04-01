process TBPROFILER_DB_UPDATE {

        /*
                This module performs TB-Profiler using TBDB results to get MT lineage and 
                identify any potential contamination in the genome
        */

        tag "$runID"

        conda params.tbprofiler_env

        container { 
                if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
                else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
                else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
                else return null
        }
        
        publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'copy'

        input:
                val(runID)

        output:
                path("update_handle.txt"), emit: tbprofiler_update_handover

        script:

                """
                # update the TBDB database
                        tb-profiler update_tbdb 
                # update the WHO database
                        tb-profiler update_tbdb --branch who

                touch update_handle.txt 
                """
}