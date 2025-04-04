process TAXONKIT_DB_UPDATE {

/*
        @author:  Poppy J Hesketh Best
        @date:    2025-04-01
        @version: 0.1
        @description: 
                This module updates the TBDB database and the WHO database using tb-profiler.
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
                path("update_db.txt"), emit: taxonkit_update_db
                path("taxonkit/*")

        script:

        """
        mkdir -p taxonkit/

        # update the TBDB database
        wget http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -xzf taxdump.tar.gz

        # move to output direcotort
        mv *.dmp taxonkit/
        mv *.prt taxonkit/
        mv *.txt taxonkit/
        rm taxdump.tar.gz

        touch update_db.txt
        """
}

/*
        >> update_db.txt 2>&1
*/