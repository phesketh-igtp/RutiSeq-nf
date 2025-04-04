process INSPECT_BBDD {

/*
    @author:  Poppy J Hesketh Best
    @date:    2025-04-01
    @version: 0.1
    @description: 
        This process checks the BBDD for the presence of negative controls outputs and generates a controls.tuple.csv file that
        is then used to create a channel of tuples with the sampleID and the paths to the forward and reverse reads. If the
        channel has no reads, it will be empty and the sampleID will not be pased onto the modules for analysis (MTBseq+TBProfiler).
        The process will not fail if the files do not exist.
*/

    tag "${sampleID}"

    array 50

    input:
        tuple val(sampleID), 
            path(forward), 
            path(reverse)

    output:
        path("controls.tuple.csv"), emit: controls_paths

    script:
        def forward_path    = forward.toRealPath()
        def reverse_path    = reverse.toRealPath()
        def cn_tbprof       = "${params.outdir}/negative-controls/tbprofiler/${sampleID}_tb_profiler.log"
        def cn_tbprof_log   = "${params.outdir}/negative-controls/tbprofiler/results/tbdb-${sampleID}.results.txt"
        def cn_mtbseq       = "${params.outdir}/negative-controls/mtbseq/Statistics/${sampleID}.Strain_Classification.tab"
        def cn_taxonomy     = "${params.outdir}/negative-controls/results/Classification/${sampleID}.k2.report"
        def cn_stats        = "${params.outdir}/negative-controls/results/Statistics/${sampleID}.stats.tsv"

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
