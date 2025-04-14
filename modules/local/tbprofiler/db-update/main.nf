process TBPROFILER_DB_UPDATE {

    /*
        @author:  Poppy J Hesketh Best
        @date:    2025-04-01
        @version: 1.2.0
        @description: 
                This module updates the TBDB database and the WHO database using tb-profiler if 
                the last update was more than a week ago.
        @changelog:
                v1.0.0-2024-04-01: Initial version
                v1.1.0-2024-04-01: Created 'log/' directory to ouput the ${sampleID}_tbprofiler.log
                v1.2.0-2025-04-01: Added check to only update if last update was more than a week ago
    */

        tag "$runID"

        conda params.tbprofiler_env

        container { 
                if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
                else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
                else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
                else return null
        }

        publishDir "${params.outDir}/db/tbprofiler/", mode: 'copy', overwrite: true

        input:
                val(runID)

        output:
                path("update_db.txt"), emit: tbprofiler_update_db
                path("tbdb/*"), optional: true

        script:

        """
                # Check if last_update.txt exists
                        if [ -f ${params.outDir}/db/tbprofiler/last_update.txt ]; then
                                last_update=\$(cat ${params.outDir}/db/tbprofiler/last_update.txt)
                                current_time=\$(date +%s)
                                week_in_seconds=604800
                        
                # Calculate the difference between current time and last update
                                time_diff=\$((current_time - last_update))

                # If the difference is less than a week, exit
                                if [ \$time_diff -lt \$week_in_seconds ]; then
                                        echo "Last update was less than a week ago. Skipping update."
                                        tb-profiler update_tbdb
                                        tb-profiler update_tbdb --branch who
                                        touch update_db.txt
                                        exit 0
                                fi
                        fi

                # If we reach here, we need to update the databases
                tb-profiler update_tbdb
                tb-profiler update_tbdb --branch who

                # update the TBDB database
                        tb-profiler update_tbdb \\
                                --branch tbdb \\
                                --db_dir .
                                
                # update the WHO database
                        tb-profiler update_tbdb \\
                                --branch who \\
                                --db_dir .

                # Record the current time as the last update time
                        date +%s > ${params.outDir}/db/tbprofiler/last_update.txt

                        touch update_db.txt
        """
}