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
        val(sampleID_list)
        //path(tbprofiler_results)

    output:
    tuple path("tbdb-tbprofiler.txt"), 
        path("who-tbprofiler.txt"),     
        path("lineages.fractions.txt"), emit: tbdb_out

    script:
        """
        mkdir -p results/; mkdir -p bam/; mkdir -p vcf/; mkdir -p tmp/

        # create the symbolic links to the result directories
            ln -s ${params.outDir}/db/samples/*/tbprofiler/tbdb-* results/
            tb-profiler collate
            sed 's/tbdb-//g' tbprofiler.txt > tbdb-tbprofiler.txt

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
            ln -s ${params.outDir}/db/samples/*/tbprofiler/who-* results/
            tb-profiler collate
            sed 's/who-//g' tbprofiler.txt > who-tbprofiler.txt
        """
}