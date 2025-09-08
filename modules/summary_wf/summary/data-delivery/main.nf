process DATA_DELIVERY{

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
            if [[ -f "${params.outDir}/${params.outDir}/negative-controls/negative-controls.xlsx" ]]; then
                cp ${params.outDir}/${params.outDir}/negative-controls/negative-controls.xlsx ${params.outDir}/results/${runID}/controls/
            fi

        # TBProfiler results 
            mkdir -p ${params.outDir}/results/${runID}/tbprofiler/
            if [[ -f "${params.outDir}/bbdd/tbprofiler/tbdb-tbprofiler.txt" ]]; then
                cp ${params.outDir}/bbdd/tbprofiler/tbdb-tbprofiler.txt ${params.outDir}/results/${runID}/tbprofiler/tbdb-tbprofiler.txt
            fi
            if [[ -f "${params.outDir}/bbdd/tbprofiler/who-only/who-tbprofiler.txt" ]]; then
                cp ${params.outDir}/bbdd/tbprofiler/who-only/who-tbprofiler.txt ${params.outDir}/results/${runID}/tbprofiler/who-tbprofiler.txt
            fi

        # MTBseq results
            mkdir -p ${params.outDir}/results/${runID}/mtbseq/
            if [[ -f "${params.outDir}/bbdd/mtbseq/pairwise/Mapping_and_Variant_Statistics.tab" ]]; then
                cp ${params.outDir}/bbdd/mtbseq/pairwise/Mapping_and_Variant_Statistics.tab ${params.outDir}/results/${runID}/mtbseq/
            fi
            if [[ -f "${params.outDir}/bbdd/mtbseq/pairwise/Strain_Classification.tab" ]]; then
                cp ${params.outDir}/bbdd/mtbseq/pairwise/Strain_Classification.tab ${params.outDir}/results/${runID}/mtbseq/
            fi

        # Matrices
            mkdir -p ${params.outDir}/results/${runID}/matrices/
            if compgen -G "${params.outDir}/bbdd/mtbseq/pairwise/*/Matrices/*.matrix.tsv.gz" > /dev/null; then
                cp ${params.outDir}/bbdd/mtbseq/pairwise/*/Matrices/*.matrix.tsv ${params.outDir}/results/${runID}/matrices/
            fi

        # Tidy up the phylogeny output
            mkdir -p ${params.outDir}/results/${runID}/phylogeny/html-out
            mkdir -p ${params.outDir}/results/${runID}/phylogeny/Rdata-out

        # Move PDF files if they exist
            if compgen -G "${params.outDir}/results/${runID}/phylogeny/*.html" > /dev/null; then
                mv ${params.outDir}/results/${runID}/phylogeny/*.pdf ${params.outDir}/results/${runID}/phylogeny/html-out/
            fi

        # Move Rdata files if they exist
            if compgen -G "${params.outDir}/results/${runID}/phylogeny/*.Rdata" > /dev/null; then
                mv ${params.outDir}/results/${runID}/phylogeny/*.Rdata ${params.outDir}/results/${runID}/phylogeny/Rdata-out/
            fi

        # Clean up: remove litter from the nexus generation
            if compgen -G "${params.outDir}/results/${runID}/networks/*join*tab" > /dev/null; then
                rm -rf ${params.outDir}/results/${runID}/networks/*join*tab
            fi
            if [[ -f "${params.outDir}/results/${runID}/phylogeny/nexus.TT.tuple.csv" ]]; then
                rm -f ${params.outDir}/results/${runID}/phylogeny/nexus.TT.tuple.csv
            fi
            if [[ -d "${params.outDir}/results/${runID}/snps/cleanup-handover" ]]; then
                rm -rf ${params.outDir}/results/${runID}/snps/cleanup-handover
            fi
        """

}