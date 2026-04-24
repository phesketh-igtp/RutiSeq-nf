process SNIPPY_ENV {
    conda params.tbprofiler_env

    container { 
        if (workflow.containerEngine == 'singularity') {
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
        } else { 
            'quay.io/biocontainers/tb-profiler' 
        }
    }

    storeDir "${params.outdir}/version_control"
    storeDir "${params.outDir}/results/${runID}/version_control/"

    input:
        val runID

    output:
        file("snippy_env_${runID}.yml")

    script:

    """
    conda export > snippy_env_${runID}.yml
    """
}