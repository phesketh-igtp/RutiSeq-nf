process FILE_CHECK {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.1.0
    @description: 
        This process checks for the existence of specific files related to MTBseq and TBProfiler results.
        It generates a CSV file with the sample ID and paths to the forward and reverse reads, as well 
        as the paths to the MTBseq and TBProfiler results if they exist.
    @changelog:
        2025-04-01: Addition of comments
        2025-11-12: Restructuring DB format to be sample-first, and results nested by samples
        2025-11-14: Improved error handling, better shell logic, fixed array syntax
*/

    tag "$sampleID"

    maxForks 200  // Changed from 'array' (not a valid directive) to 'maxForks'

    input:
        tuple val(sampleID), 
            path(forward), 
            path(reverse), 
            val(type)

    output:
        path("${sampleID}.tuple.csv"), emit: sample_paths

    script:
        // Handle empty forward/reverse (for complete samples)
        def forward_path = forward ? forward.toRealPath() : 'null'
        def reverse_path = reverse ? reverse.toRealPath() : 'null'
        
        // Define expected result file paths
        def mtbseq_class = "${params.outDir}/db/samples/${sampleID}/mtbseq/Classification/${sampleID}.Strain_Classification.tab"
        def mtbseq_stats = "${params.outDir}/db/samples/${sampleID}/mtbseq/Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab"
        def mtbseq_pos   = "${params.outDir}/db/samples/${sampleID}/mtbseq/Position_Tables/${sampleID}.gatk_position_table.tab"
        def mtbseq_vars  = "${params.outDir}/db/samples/${sampleID}/mtbseq/Called/${sampleID}.gatk_position_variants_cf4_cr4_fr75_ph4_outmode000.tab"
        def tbdb_out     = "${params.outDir}/db/samples/${sampleID}/tbprofiler/tbdb-${sampleID}.results.txt"
        def who_out      = "${params.outDir}/db/samples/${sampleID}/tbprofiler/who-${sampleID}.results.txt"
        def snippy_vcf   = "${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.vcf"

        """
        #!/bin/bash
        set -euo pipefail

        # Check if ALL result files exist
        if [[ -f "${mtbseq_class}" ]] && \
            [[ -f "${mtbseq_stats}" ]] && \
            [[ -f "${mtbseq_pos}" ]] && \
            [[ -f "${mtbseq_vars}" ]] && \
            [[ -f "${tbdb_out}" ]] && \
            [[ -f "${who_out}" ]] && \
            [[ -f "${snippy_vcf}" ]]; then
            
            # Sample is complete - use result files, set FASTQs to null
            echo "${sampleID},null,null,${type},${mtbseq_class},${mtbseq_stats},${mtbseq_pos},${mtbseq_vars},${tbdb_out},${who_out},${snippy_vcf}" > "${sampleID}.tuple.csv"
            echo "[FILE_CHECK] ${sampleID}: Complete (using existing results)" >&2
        else
            # Sample needs processing - use FASTQ paths, set results to null
            echo "${sampleID},${forward_path},${reverse_path},${type},null,null,null,null,null,null,null" > "${sampleID}.tuple.csv"
            echo "[FILE_CHECK] ${sampleID}: Incomplete (will process FASTQs)" >&2
        fi
        """
}