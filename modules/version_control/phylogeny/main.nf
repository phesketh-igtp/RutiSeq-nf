process PHYLOGENY_ENV {

    conda params.phylogeny_env

    container { 
        if (workflow.containerEngine == 'singularity') {
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
        } else { 
            'quay.io/biocontainers/tb-profiler' 
        }
    }

    storeDir "${params.outDir}/results/${params.runID}/version_control/"

    input:
        val runID

    output:
        file("phylogeny_env_${runID}.yml")

    script:

    """
    conda export > phylogeny_env_${runID}.yml
    """
}