process PREPARE_DATA_DELIVERY{

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0.0
    @description:
        This process prepares the data delivery for the final results. It moves the relevant files to the appropriate directories and cleans up any unnecessary files.
        It is called after the summary report generation and before the final data delivery.
    @changelog
        v1.0.0-2025-04-01: Initial version
*/

    tag "${runID}"

    input:
        val(runID)
        path(handover)

    script:
        """
        # Controls results
        mkdir -p ${params.outDir}/results/${runID}/controls/
            cp ${params.outDir}/${params.outDir}/negative-controls/negative-controls.xlsx ${params.outDir}/results/${runID}/controls/

        # TBProfiler results 
            mkdir -p ${params.outDir}/results/${runID}/tbprofiler/
            cp ${params.outDir}/bbdd/tbprofiler/tbdb-tbprofiler.txt ${params.outDir}/results/${runID}/tbprofiler/tbdb-tbprofiler.txt
            cp ${params.outDir}/bbdd/tbprofiler/who-only/who-tbprofiler.txt ${params.outDir}/results/${runID}/tbprofiler/who-tbprofiler.txt

        # MTBseq results
        mkdir -p ${params.outDir}/results/${runID}/mtbseq/
            cp ${params.outDir}/bbdd/mtbseq/pairwise/Mapping_and_Variant_Statistics.tab ${params.outDir}/results/${runID}/mtbseq/Mapping_and_Variant_Statistics.tab
            cp ${params.outDir}/bbdd/mtbseq/pairwise/Strain_Classification.tab ${params.outDir}/results/${runID}/mtbseq/Strain_Classification.tab

        # Matrices
            mkdir -p ${params.outDir}/results/${runID}/matrices/
            cp ${params.outDir}/bbdd/mtbseq/pairwise/*/Matrices/*.matrix.tsv.gz ${params.outDir}/results/${runID}/matrices/

        # Tidy up the phylogeny output
            mkdir -p ${params.outDir}/results/${runID}/phylogeny/pdf-out
            mkdir -p ${params.outDir}/results/${runID}/phylogeny/Rdata-out

        # Move PDF files if they exist
            if ls ${params.outDir}/results/${runID}/phylogeny/*.pdf 1> /dev/null 2>&1; then
                mv ${params.outDir}/results/${runID}/phylogeny/*.pdf ${params.outDir}/results/${runID}/phylogeny/pdf-out/
            fi

        # Move Rdata files if they exist
            if ls ${params.outDir}/results/${runID}/phylogeny/*Rdata 1> /dev/null 2>&1; then
                mv ${params.outDir}/results/${runID}/phylogeny/*.Rdata ${params.outDir}/results/${runID}/phylogeny/Rdata-out/
            fi

        # Clean up: remove litter from the nexus generation
            rm -rf ${params.outDir}/results/${runID}/networks/*join*tab
            rm -rf ${params.outDir}/results/${runID}/phylogeny/nexus.TT.tuple.csv
            rm -rf ${params.outDir}/results/${runID}/snps/cleanup-handover
        """

}