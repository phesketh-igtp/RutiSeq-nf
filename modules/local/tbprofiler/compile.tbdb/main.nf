process TBPROFILER_COMPILE_TBDB {

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
        path("tbdb-tbprofiler.txt"),                 emit: tbdb_results
        path("lineages.fractions.txt"),              emit: lineage_fractions
        path("tbdb-tbprofiler.dr.indiv.itol.txt")
        path("tbdb-tbprofiler.dr.itol.txt")
        path("tbdb-tbprofiler.lineage.itol.txt")
        path("tbdb-tbprofiler.variants.csv")
        path("tbdb-tbprofiler.variants.txt")

    script:
        """
        mkdir -p results/; mkdir -p bam/; mkdir -p vcf/; mkdir -p tmp/

        # create the symbolic links to the result directories
                mv tbdb-* tmp/
            ln -s ${params.outdir}/bbdd/tbprofiler/results/* results/;
            ln -s ${params.outdir}/bbdd/tbprofiler/bam/* bam/;
            ln -s ${params.outdir}/bbdd/tbprofiler/vcf/* vcf/

        # perform tb profiler compile:
            tb-profiler collate --full --mark_missing --all_variants --itol

            sed -i 's/tbdb-//g' tbprofiler.txt
            sed -i 's/tbdb-//g' tbprofiler.variants.csv
            sed -i 's/tbdb-//g' tbprofiler.variants.txt
            sed -i 's/tbdb-//g' tbprofiler.dr.indiv.itol.txt
            sed -i 's/tbdb-//g' tbprofiler.dr.itol.txt
            sed -i 's/tbdb-//g' tbprofiler.lineage.itol.txt
            
            mv tmp/* .

        # Move the files to give them unique names
            mv tbprofiler.txt tbdb-tbprofiler.txt
            mv tbprofiler.variants.csv          tbdb-tbprofiler.variants.csv
            mv tbprofiler.variants.txt          tbdb-tbprofiler.variants.txt
            mv tbprofiler.dr.indiv.itol.txt     tbdb-tbprofiler.dr.indiv.itol.txt
            mv tbprofiler.dr.itol.txt           tbdb-tbprofiler.dr.itol.txt
            mv tbprofiler.lineage.itol.txt      tbdb-tbprofiler.lineage.itol.txt

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
        """
}