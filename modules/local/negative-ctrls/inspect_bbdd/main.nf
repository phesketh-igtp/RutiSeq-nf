process INSPECT_BBDD {

/*
    @author:  Poppy J Hesketh Best
    @date:    2025-04-07
    @version: 1.0.0
    @description: 
        This process checks the BBDD for the presence of negative controls outputs and generates a controls.tuple.csv file that
        is then used to create a channel of tuples with the sampleID and the paths to the forward and reverse reads. If the
        channel has no reads, it will be empty and the sampleID will not be pased onto the modules for analysis (MTBseq+TBProfiler).
        The process will not fail if the files do not exist.
    @changelog:
        v1.0.0-2025-04-01:  Added - Initial version
        v1.0.1-2025-04-07:  Fixed - Corrected path for tbprofiler log, k2 taxonomy, and stats
*/

    tag "${sampleID}"

    array 100

    input:
        tuple val(sampleID), 
            path(forward), 
            path(reverse)

    output:
        path("controls.tuple.csv"), emit: controls_paths

    script:
        def forward_path    = forward.toRealPath()
        def reverse_path    = reverse.toRealPath()
        def cn_tbprof       = "${params.outDir}/negative-controls/tbprofiler/log/${sampleID}_tb_profiler_status.txt"
        def cn_tbprof_log   = "${params.outDir}/negative-controls/tbprofiler/results/tbdb-${sampleID}.results.txt"
        def cn_mtbseq       = "${params.outDir}/negative-controls/mtbseq/Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab"
        def cn_taxonomy     = "${params.outDir}/negative-controls/Classification/${sampleID}.k2.report"
        def cn_stats        = "${params.outDir}/negative-controls/Statistics/${sampleID}.stats.tsv"

        """
        # Check if the files exist and create the controls.tuple.csv file
        if ([[ -f "${cn_tbprof}" ]] || [[ -f "${cn_tbprof_log}" ]]) && [[ -f "${cn_mtbseq}" ]] && [[ -f "${cn_taxonomy}" ]] && [[ -f "${cn_stats}" ]]; then
            echo "${sampleID},," > controls.tuple.csv
            echo "DEBUG: Added to samples.txt: ${sampleID}" >&2
        else
            echo "${sampleID},${forward_path},${reverse_path}" > controls.tuple.csv
            echo "DEBUG: Added to samples.txt: ${sampleID}" >&2
        fi

        # Ensure the process doesn't fail if one of the files doesn't exist
            touch controls.tuple.csv
        """
}
