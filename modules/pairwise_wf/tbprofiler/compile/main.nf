process TBPROFILER_COMPILE {

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

    container { 
        if (workflow.containerEngine == 'singularity') {
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
        } else { 
            'quay.io/biocontainers/tb-profiler' 
        }
    }

    publishDir "${params.outDir}/db/comparison/tbprofiler/", mode: 'copy'

    input:
        path(tbprofiler_results)

    output:
        path("tbdb-tbprofiler.txt"),                 emit: tbdb_results
        path("lineages.fractions.txt"),              emit: lineage_fractions
        path("tbdb-tbprofiler.dr.indiv.itol.txt")
        path("tbdb-tbprofiler.dr.itol.txt")
        path("tbdb-tbprofiler.lineage.itol.txt")
        path("tbdb-tbprofiler.variants.csv")
        path("tbdb-tbprofiler.variants.txt")
        path("who-tbprofiler.txt"),                 emit: who_results
        path("who-tbprofiler.dr.indiv.itol.txt")
        path("who-tbprofiler.dr.itol.txt")
        path("who-tbprofiler.lineage.itol.txt")
        path("who-tbprofiler.variants.csv")
        path("who-tbprofiler.variants.txt")

    script:
        """
        mkdir -p results/; mkdir -p bam/; mkdir -p vcf/; mkdir -p tmp/

        # create the symbolic links to the result directories
            mv tbdb-* tmp/
            ln -s ${params.outDir}/db/*/tbprofiler/results/tbdb-* results/

        # perform tb profiler compile:
            tb-profiler collate --full --mark_missing --all_variants --itol

            sed 's/tbdb-//g' tbprofiler.txt                 > tbdb-tbprofiler.txt
            sed 's/tbdb-//g' tbprofiler.variants.csv        > tbdb-tbprofiler.variants.csv
            sed 's/tbdb-//g' tbprofiler.variants.txt        > tbdb-tbprofiler.variants.txt
            sed 's/tbdb-//g' tbprofiler.dr.indiv.itol.txt   > tbdb-tbprofiler.dr.indiv.itol.txt
            sed 's/tbdb-//g' tbprofiler.dr.itol.txt         > tbdb-tbprofiler.dr.itol.txt
            sed 's/tbdb-//g' tbprofiler.lineage.itol.txt    > tbdb-tbprofiler.lineage.itol.txt

        # Get the fractions of all the lineages
            for file in results/tbdb-*.txt; do
                id=\$(basename \$file .results.txt)
                sed -e '/Resistance report/,\$d' \\
                    -e '1,/Lineage report/d' \\
                    -e 's@-@@g' \\
                    -e '/^\$/d' \\
                    -e "s@^@\${id}\t@" \\
                    \$file | sort >> lineages.fractions.txt
            done

        # Add header to lineages.fractions.txt
            sed -i '1iSampleID\\tLineage\\tFraction\\tFamily\\tRd' lineages.fractions.txt
            sed -i 's/tbdb-//g' lineages.fractions.txt

        # WHO database one:
            rm results/*
            ln -s ${params.outDir}/db/*/tbprofiler/results/who-* results/;

        # perform tb profiler compile:
            tb-profiler collate --full --mark_missing --all_variants --itol
            mv tbprofiler.txt tbdb-tbprofiler.txt
            mv tbprofiler.variants.csv          tbdb-tbprofiler.variants.csv
            mv tbprofiler.variants.txt          tbdb-tbprofiler.variants.txt
            mv tbprofiler.dr.indiv.itol.txt     tbdb-tbprofiler.dr.indiv.itol.txt
            mv tbprofiler.dr.itol.txt           tbdb-tbprofiler.dr.itol.txt
            mv tbprofiler.lineage.itol.txt      tbdb-tbprofiler.lineage.itol.txt
        """
}