process TBPROFILER_COMPILE_TBDB {

    tag "${runID}"

    conda params.tbprofiler_env

    container { 
        if (workflow.containerEngine == 'singularity') {
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
        } else { 
            'quay.io/biocontainers/tb-profiler' 
        }
    }

    publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'copy'

    input:
        val(runID)
        path(tbprofiler_results)

    output:
        path("tbprofiler.txt"),                 emit: tbdb_results
        path("tbprofiler.dr.indiv.itol.txt")
        path("tbprofiler.dr.itol.txt")
        path("tbprofiler.lineage.itol.txt")
        path("tbprofiler.variants.csv")
        path("tbprofiler.variants.txt")

    script:
        """
        mkdir -p results/; mkdir -p  bam/; mkdir -p vcf/

        # create the symbolic links to the result directories
        ln -s ${params.outdir}/bbdd/tbprofiler/results/* results/
        ln -s ${params.outdir}/bbdd/tbprofiler/bam/* bam/
        ln -s ${params.outdir}/bbdd/tbprofiler/vcf/* vcf/

        tb-profiler collate --full --mark_missing --all_variants --itol

        sed -i 's/tbdb-//g' tbprofiler.txt
        sed -i 's/tbdb-//g' tbprofiler.dr.indiv.itol.txt
        sed -i 's/tbdb-//g' tbprofiler.dr.itol.txt
        sed -i 's/tbdb-//g' tbprofiler.lineage.itol.txt
        sed -i 's/tbdb-//g' tbprofiler.variants.csv
        sed -i 's/tbdb-//g' tbprofiler.variants.txt
        
        """
}