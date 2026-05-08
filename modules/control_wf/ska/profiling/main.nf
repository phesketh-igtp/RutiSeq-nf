process SKA_PROFILING {

    conda params.readQC_env

    publishDir "${params.outDir}/db/qc/${params.runID}/", 
        mode: 'copy',
        overwrite: true

    input:
        path(samplesheet, stageAs: "samplesheet.csv")

    output:
        path('ska_distances.txt'), emit: ska_out
    
    script:
    
    """
    # Make the sample sheet for ska2
        sed 's@,@\t@g' samplesheet.csv \
            | cut -f2,3,4 \
            | tail -n +2 \
            | sed '/^[[:space:]]*\$/d' \
            > input.txt

    # Run SKA profiling
        ska build \\
            -o seqs -f input.txt \\
            -k ${params.ska_kmer} \\
            --min-count ${params.ska_min_count} \\
            --proportion-reads ${params.ska_proportion_reads} \\
            --qual-filter ${params.ska_qual_filter} \\
            --min-qual ${params.ska_min_qual} \\
            --threads ${task.cpus} \\
            --verbose \\

        ska distance ${params.ska_distance_args} \\
            -o ska_distances.txt \\
            seqs.skf \\
            --min-freq ${params.ska_min_freq}
    """
}


/*
@author: Poppy J hesketh Best
@date: 2025-11-17
@version: v1.2.0
@description:
    Created an input file from the samplesheet that is used by ska2
        to generate the ska profiles, before creating a 
        distance matrix of those profiles.
@changelog:
    v1.0.0-2025-11-17: Functioning module created.
    v1.1.0-2025-11-26: Modified to add params for ska2 profiling (--verbose).
    v1.2.0-2025-12-19: Modified to handle when the samlesheet has an empty line
*/