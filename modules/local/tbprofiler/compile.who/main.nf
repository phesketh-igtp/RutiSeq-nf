process TBPROFILER_COMPILE_WHO {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-11
    @version: 1.1.0
    @description:
        This process compiles the TB-Profiler results from the tbdb pipeline
        into a single file. Since TB-Profiler requires the results directory to be 
        present in the current working directory, we create symbolic links to the
        results directory, the bam directory and the vcf directory. The symbolic
        links are then used to run the tb-profiler collate command. The results
        are then moved to the current working directory and renamed to remove the
        tbdb- prefix. Renaming is to prevent clashes with input files in downstream
        processes. The results are then moved to the output directory.
    @changelog:
        v1.0.0-2024-12-01: Initial version added
        v1.1.0-2025-04-11: Added - a handover from the TBPROFILER db updaitng module
*/

    conda params.tbprofiler_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
            } else { 'quay.io/biocontainers/tb-profiler' }
    }
    
    publishDir "${params.outdir}/bbdd/tbprofiler/who-only/", mode: 'copy'

    input:
        val runID
        path (tbprofiler_who_results)
        path(tbprofiler_update_db)

    output:
        path("who-tbprofiler.txt"),               emit: who_results
        path("who-tbprofiler.variants.csv")
        path("who-tbprofiler.variants.txt")

    script:

        """
        mkdir -p results/; mkdir -p bam/; mkdir -p vcf/

        # create the symbolic links to the result directories
        ln -s ${params.outdir}/bbdd/tbprofiler/who-only/results/* results/
        #ln -s ${params.outdir}/bbdd/tbprofiler/who-only/bam/* bam/
        #ln -s ${params.outdir}/bbdd/tbprofiler/who-only/vcf/* vcf/

        tb-profiler collate # --full --mark_missing --all_variants --itol

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