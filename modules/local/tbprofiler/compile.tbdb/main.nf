process TBPROFILER_COMPILE_TBDB {
    tag "${runID}"

    conda 'bioconda::tb-profiler==6.5.0'

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
        path("sample_lineage.txt"),             emit: sample_lineage
        path("tbprofiler.txt")
        path("tbprofiler.dr.indiv.itol.txt")
        path("tbprofiler.dr.itol.txt")
        path("tbprofiler.lineage.itol.txt")
        path("tbprofiler.variants.csv")
        path("tbprofiler.variants.txt")

    script:
        """
        mkdir -p results/
        mkdir -p bam/
        mkdir -p vcf/

        tb-profiler collate --full --mark_missing --all_variants --itol

        cut -f1,3 tbprofiler.txt | sed '1d' > tuple_lineages.tsv

        # Extract sample ID and lineage
        SAMPLE_LINEAGE=\$(head -n 1 tuple_lineages.tsv)
        IFS=\$'\\t' read -r SAMPLE_ID LINEAGE <<< "\$SAMPLE_LINEAGE"
        
        if [[ -z "\$SAMPLE_ID" ]]; then
            SAMPLE_ID="NO_SAMPLE"
        fi
        
        if [[ -z "\$LINEAGE" ]]; then
            LINEAGE="0"
        fi

        # Output the sample ID and lineage for Nextflow to capture
        echo "\$SAMPLE_ID\\t\$LINEAGE" > sample_lineage.txt
        """
}