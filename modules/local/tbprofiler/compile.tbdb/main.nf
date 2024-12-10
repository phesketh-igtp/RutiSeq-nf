process TBPROFILER_COMPILE_TBDB {
    
    tag "${runID}"

    conda 'bioconda::tb-profiler==6.5.0'

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
            } else { 'quay.io/biocontainers/tb-profiler' }
    }

    publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'copy'

    input:
        val runID

    output:
    path("tbprofiler.txt"),             emit: tbprofile_who_compile
    path "tbprofiler.dr.indiv.itol.txt"
    path "tbprofiler.dr.itol.txt"
    path "tbprofiler.lineage.itol.txt"
    path "tbprofiler.variants.csv"
    path "tbprofiler.variants.txt"

    script:

    """
    mkdir results/

    ln -s ${params.outdir}/bbdd/tbprofiler/results/* results/

    tb-profiler collate --full --mark_missing --all_variants --itol

    rm -rf results/

    """
    
}