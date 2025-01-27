process INSPECT_BBDD {

    tag "$sampleID"

    array 10

    input:
        tuple val(sampleID), path(forward), path(reverse)

    output:
        path("controls.tuple.csv"), emit: controls_paths

    script:
        def forward_path = forward.toRealPath()
        def reverse_path = reverse.toRealPath()
        def qc_results   = "${params.outdir}/bbdd/negative-controls/${sampleID}.qc.out"

        """
        if [ -f ${qc_results} ]; then       
            echo "${sampleID},,,${qc_results}" > controls.tuple.csv
            echo "DEBUG: Added to samples.txt: ${sampleID}" >&2    
        else
            echo "${sampleID},${forward_path},${reverse_path}," > controls.tuple.csv
            echo "DEBUG: Added to samples.txt: ${sampleID}" >&2
        fi

        # Ensure the process doesn't fail if one of the files doesn't exist
        touch controls.tuple.csv

        """
}
