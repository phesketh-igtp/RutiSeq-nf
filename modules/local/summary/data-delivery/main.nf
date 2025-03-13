process PREPARE_DATA_DELIVERY{

    tag "${runID}"

    input:
        val(runID)
        path(handover)

    script:
        """
        # TBProfiler results 
        mkdir -p ${params.outdir}/results/${runID}/tbprofiler/
        cp ${params.outdir}/bbdd/tbprofiler/tbdb-tbprofiler.txt ${params.outdir}/results/${runID}/tbprofiler/tbdb-tbprofiler.txt
        cp ${params.outdir}/bbdd/tbprofiler/who-only/who-tbprofiler.txt ${params.outdir}/results/${runID}/tbprofiler/who-tbprofiler.txt

        # MTBseq results
        mkdir -p ${params.outdir}/results/${runID}/mtbseq/
        cp ${params.outdir}/bbdd/mtbseq/pairwise/Mapping_and_Variant_Statistics.tab ${params.outdir}/results/${runID}/mtbseq/Mapping_and_Variant_Statistics.tab
        cp ${params.outdir}/bbdd/mtbseq/pairwise/Strain_Classification.tab ${params.outdir}/results/${runID}/mtbseq/Strain_Classification.tab

        # Matrices
        mkdir -p ${params.outdir}/results/${runID}/matrices/
        cp ${params.outdir}/bbdd/mtbseq/pairwise/*/Matrices/*.matrix.tsv.gz ${params.outdir}/results/${runID}/matrices/

        # Tidy up the phylogeny output
        mkdir -p ${params.outdir}/results/${runID}/phylogeny/pdf-out
        mkdir -p ${params.outdir}/results/${runID}/phylogeny/Rdata-out
        mv ${params.outdir}/results/${runID}/phylogeny/*.pdf ${params.outdir}/results/${runID}/phylogeny/pdf-out/
        mv ${params.outdir}/results/${runID}/phylogeny/*RData ${params.outdir}/results/${runID}/phylogeny/Rdata-out/

        # Clean up: remove litter from the nexus generation
        rm -rf ${params.outdir}/results/${runID}/networks/*join*tab
        rm -rf ${params.outdir}/results/${runID}/phylogeny/nexus.TT.tuple.csv
        rm -rf ${params.outdir}/results/${runID}/snps/cleanup-handover
        """

}