process TBPROFILER_DB_UPDATE {

/*
        @author:  Poppy J Hesketh Best
        @date:    2025-04-01
        @version: 0.1
        @description: 
                This module updates the TBDB database and the WHO database using tb-profiler.
*/

        tag "$runID"

        conda params.tbprofiler_env

        container { 
                if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
                else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
                else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
                else return null
        }
        
        publishDir "${params.outdir}/db/tbprofiler/", mode: 'copy', overwrite: true

        input:
                val(runID)

        output:
                path("update_db.txt"), emit: tbprofiler_update_db
                path("tbdb/*")

        script:

        """
        mkdir -p ${params.outdir}/db/tbprofiler/ 

        # update the TBDB database
                tb-profiler update_tbdb \
                --db_dir . \
                > update_db.txt 2>&1

        # update the WHO database
                tb-profiler update_tbdb \
                --branch who \
                --db_dir . \
                >> update_db.txt 2>&1
        """
}