process INSPECT_BBDD {

    tag "$sampleID"

    array 10

    input:
        tuple val(sampleID), 
            path(forward), 
            path(reverse)

    output:
        path("controls.tuple.csv"), emit: controls_paths

    script:
        def forward_path    = forward.toRealPath()
        def reverse_path    = reverse.toRealPath()
        def cn_tbprof       = "${params.outdir}/bbdd/negative-controls/tbprofiler/${sampleID}_tb_profiler.log"
        def cn_tbprof_log   = "${params.outdir}/bbdd/negative-controls/tbprofiler/results/tbdb-${sampleID}.results.txt"
        def cn_mtbseq       = "${params.outdir}/bbdd/negative-controls/mtbseq/${sampleID}"

        """
        if { [ -f "${cn_tbprof}" ] || [ -f "${cn_tbprof_log}" ]; } && [ -f "${cn_mtbseq}" ]; then
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
