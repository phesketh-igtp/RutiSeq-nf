process DATA_DELIVERY {

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

    tag "${params.runID}"

    input:
        path(sylph_results)

        tuple path(html_report, stageAs: "${params.runID}_reads-controls-report.html"),
            path(warnings, stageAs: 'warnings.out')
        
        path(nexus_handover)
        

    script:

        """
        # Report output:
            cp ${params.runID}_reads-controls-report.html ${params.outDir}/results/${params.runID}/${params.runID}_reads-controls-report.html

        # TBProfiler results 
            mkdir -p ${params.outDir}/results/${params.runID}/tbprofiler/
            if [[ -f "${params.outDir}/db/comparison/tbprofiler/tbdb-tbprofiler.txt" ]]; then
                cp ${params.outDir}/db/comparison/tbprofiler/tbdb-tbprofiler.txt ${params.outDir}/results/${params.runID}/tbprofiler/tbdb-tbprofiler.txt
            fi
            if [[ -f "${params.outDir}/db/comparison/tbprofiler/who-tbprofiler.txt" ]]; then
                cp ${params.outDir}/db/comparison/tbprofiler/who-tbprofiler.txt ${params.outDir}/results/${params.runID}/tbprofiler/who-tbprofiler.txt
            fi

        # MTBseq results
            mkdir -p ${params.outDir}/results/${params.runID}/mtbseq/
            if [[ -f "${params.outDir}/db/comparison/mtbseq/Mapping_and_Variant_Statistics.tab" ]]; then
                cp ${params.outDir}/db/comparison/mtbseq/Mapping_and_Variant_Statistics.tab ${params.outDir}/results/${params.runID}/mtbseq/
            fi
            if [[ -f "${params.outDir}/db/comparison/mtbseq/Strain_Classification.tab" ]]; then
                cp ${params.outDir}/db/comparison/mtbseq/Strain_Classification.tab ${params.outDir}/results/${params.runID}/mtbseq/
            fi

        # Matrices
            mkdir -p ${params.outDir}/results/${params.runID}/matrices/
            if compgen -G "${params.outDir}/db/comparison/mtbseq/*/Matrices/*.matrix.tsv.gz" > /dev/null; then
                cp ${params.outDir}/db/comparison/mtbseq/*/Matrices/*.matrix.tsv ${params.outDir}/results/${params.runID}/matrices/
            fi

        # Clean up: remove litter from the nexus generation
            if compgen -G "${params.outDir}/results/${params.runID}/networks/*tab" > /dev/null; then
                rm -rf ${params.outDir}/results/${params.runID}/networks/*tab
            fi
            if [[ -f "${params.outDir}/results/${params.runID}/phylogeny/nexus.TT.tuple.csv" ]]; then
                rm -f ${params.outDir}/results/${params.runID}/phylogeny/nexus.TT.tuple.csv
            fi
            if [[ -d "${params.outDir}/results/${params.runID}/snps/cleanup-handover" ]]; then
                rm -rf ${params.outDir}/results/${params.runID}/snps/cleanup-handover
            fi
        """

}