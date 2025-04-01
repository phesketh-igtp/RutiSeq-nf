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
        
        publishDir "${params.outdir}/db/tbprofiler/", mode: 'copy'

        input:
                val(runID)

        output:
                path("update_db.txt"), emit: tbprofiler_update_db

        script:

        """
        mkdir -p ${params.outdir}/db/tbprofiler/ 

        # update the TBDB database
                tb-profiler update_tbdb \
                --db_dir ${params.outdir}/db/tbprofiler/ \
                > update_db.txt 2>&1

        # update the WHO database
                tb-profiler update_tbdb \
                --branch who \
                --db_dir ${params.outdir}/db/tbprofiler/ \
                >> update_db.txt 2>&1
        """
}