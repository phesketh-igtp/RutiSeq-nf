process TBPROFILER_COMPILE_WHO {
    
    tag "$params.runID"

    conda params.tbprofiler_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
            } else { 'quay.io/biocontainers/tb-profiler' }
    }
    
    publishDir "${params.outdir}/bbdd/tbprofiler/who-only/", mode: 'copy'

    input:
        val runID
        path (tbprofiler_who_results)

    output:
        path("who-tbprofiler.txt"),               emit: who_results
        path("who-tbprofiler.variants.csv")
        path("who-tbprofiler.variants.txt")
        //path("who-tbprofiler.dr.indiv.itol.txt")
        //path("who-tbprofiler.dr.itol.txt")
        //path("who-tbprofiler.lineage.itol.txt")


    script:

        """
        mkdir -p results/; mkdir -p  bam/; mkdir -p vcf/

        # create the symbolic links to the result directories
        ln -s ${params.outdir}/bbdd/tbprofiler/results/* results/
        #ln -s ${params.outdir}/bbdd/tbprofiler/bam/* bam/
        #ln -s ${params.outdir}/bbdd/tbprofiler/vcf/* vcf/

        tb-profiler collate

        # rm the prefix in the files to enable merging later
        sed -i 's/who-//g' tbprofiler.txt
        sed -i 's/who-//g' tbprofiler.variants.csv
        sed -i 's/who-//g' tbprofiler.variants.txt
        #sed -i 's/who-//g' tbprofiler.dr.indiv.itol.txt
        #sed -i 's/who-//g' tbprofiler.dr.itol.txt
        #sed -i 's/who-//g' tbprofiler.lineage.itol.txt

        # Move the files to give them unique names
        mv tbprofiler.txt                   who-tbprofiler.txt
        mv tbprofiler.variants.csv          who-tbprofiler.variants.csv
        mv tbprofiler.variants.txt          who-tbprofiler.variants.txt
        #mv tbprofiler.dr.indiv.itol.txt    who-tbprofiler.dr.indiv.itol.txt
        #mv tbprofiler.dr.itol.txt          who-tbprofiler.dr.itol.txt
        #mv tbprofiler.lineage.itol.txt     who-tbprofiler.lineage.itol.txt

        """


}