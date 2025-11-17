process SYLPH_DB {

/*
    @author: Poppy J hesketh Best
    @date: 2025-11-17
    @version: v1.0.0
    @description:
        Download Sylph database (gtdb-r220) and stores in the storeDir, meaning
            if it exists in the storeDir, then it will not download the files.
        Outputs a channel of the database that are used by the respective tools.
    @changelog:
        v1.0.0-2025-11-17: Functioning module created.
*/

    tag "Sylph database: GTDB-R220"
    
    conda params.taxonomy_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.docker_sylph
            else if (workflow.containerEngine == 'docker') return params.docker_sylph
            else if (workflow.containerEngine == 'apptainer') return params.docker_sylph
            else return null
        }

    storeDir "${params.storeDir}/sylph/"

    input:
        val(runID)

    output:
        tuple path("gtdb-r220-c200-dbv1.syldb"),
            path("sylph-tax/"), emit: db
            path("versions.yml")

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # Download sylph database
    wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb

    mkdir -p sylph-tax/
    sylph-tax download --download-to sylph-tax/

        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    runID: ${runID}
        sylph: \$(sylph -V | cut -f2 -d ' ')
        
    END_VERSIONS
    """
}