process TAXONKIT_DB_UPDATE {

/*
    @author:  Poppy J Hesketh Best
    @date:    2025-04-04
    @version: 1.1.0
    @description: 
        This module updates the NCBI taxonomy list database for
        taxonkit to use in reformating the taxonomy names by kraken2.
        It only updates if the last update was more than a week ago.
    @dependencies:
        - taxonkit=0.19.0
    @changelog:
        v1.0.0-2025-04-04: Initial version
        v1.0.1-2025-04-05: Renamed the handover file to 'update_taxonkit_db.txt'
        v1.1.0-2025-04-08: Added check to only update if last update was more than a week ago
*/

    tag "$runID"

    conda "bioconda::taxonkit=0.19.0"

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
            else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
            else return null
        }
        
    publishDir "${params.outdir}/db/taxonkit/", mode: 'copy', overwrite: true

    input:
                val(runID)

    output:
        path("update_taxonkit_db.txt"), emit: taxonkit_update_db
        path("taxonkit/*"), optional: true

    script:

        """
        # Check if ${params.outdir}/db/taxonkit/last_update.txt exists
            if [ -f ${params.outdir}/db/taxonkit/last_update.txt ]; then
                last_update=\$(cat ${params.outdir}/db/taxonkit/last_update.txt)
                current_time=\$(date +%s)
                week_in_seconds=604800

        # Calculate the difference between current time and last update
            time_diff=\$((current_time - last_update))

        # If the difference is less than a week, exit
            if [ \$time_diff -lt \$week_in_seconds ]; then
                echo "Last update was less than a week ago. Skipping update."
                touch update_taxonkit_db.txt
                exit 0
                fi
            fi

        # If we reach here, we need to update the database
			mkdir -p taxonkit/

        # update the NCBI taxonomy database
            wget http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
            tar -xzf taxdump.tar.gz

        # move to output directory
            mv *.dmp taxonkit/
            mv *.prt taxonkit/
            mv *.txt taxonkit/
            rm taxdump.tar.gz

        # Record the current time as the last update time
            date +%s > ${params.outdir}/db/taxonkit/last_update.txt

        touch update_taxonkit_db.txt
                """
}