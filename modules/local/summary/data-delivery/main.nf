process PREPARE_DATA_DELIVERY{

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process prepares the data delivery for the final results. It moves the relevant files to the appropriate directories and cleans up any unnecessary files.
        It is called after the summary report generation and before the final data delivery.
*/

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

        # Move PDF files if they exist
            if ls ${params.outdir}/results/${runID}/phylogeny/*.pdf 1> /dev/null 2>&1; then
                mv ${params.outdir}/results/${runID}/phylogeny/*.pdf ${params.outdir}/results/${runID}/phylogeny/pdf-out/
            fi

        # Move Rdata files if they exist
            if ls ${params.outdir}/results/${runID}/phylogeny/*Rdata 1> /dev/null 2>&1; then
                mv ${params.outdir}/results/${runID}/phylogeny/*.Rdata ${params.outdir}/results/${runID}/phylogeny/Rdata-out/
            fi

        # Clean up: remove litter from the nexus generation
            rm -rf ${params.outdir}/results/${runID}/networks/*join*tab
            rm -rf ${params.outdir}/results/${runID}/phylogeny/nexus.TT.tuple.csv
            rm -rf ${params.outdir}/results/${runID}/snps/cleanup-handover
        """

}