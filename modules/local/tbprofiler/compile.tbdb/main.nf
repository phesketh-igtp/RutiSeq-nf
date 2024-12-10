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
        path("tbprofiler.txt"),                             emit: tbprofile_tdb_compile
        path "tbprofiler.dr.indiv.itol.txt"
        path "tbprofiler.dr.itol.txt"
        path "tbprofiler.lineage.itol.txt"
        path "tbprofiler.variants.csv"
        path "tbprofiler.variants.txt"

    script:

    """

    mkdir -p results/; mkdir -p  bam/; mkdir -p vcf/

    for file in ${params.outdir}/bbdd/tbprofiler/results/*; do
        ln -s "\$file" results/
    done

    for file in ${params.outdir}/bbdd/tbprofiler/vcf/*; do
        ln -s "\$file" vcf/
    done

    for file in ${params.outdir}/bbdd/tbprofiler/bam/*; do
        ln -s "\$file" bam/
    done

    for file in ${params.outdir}/bbdd/tbprofiler/results/*.results.json; do
        ln -s "\$file" results/
    done

    tb-profiler collate --full --mark_missing --all_variants --itol

    """

}