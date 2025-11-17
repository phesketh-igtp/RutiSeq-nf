process SKA_PROFILING {

/*
    @author: Poppy J hesketh Best
    @date: 2025-11-17
    @version: v1.0.0
    @description:
        Created an input file from the samplesheet that is used by ska2
            to generate the ska profiles, before creating a 
            distance matrix of those profiles.
    @changelog:
        v1.0.0-2025-11-17: Functioning module created.
*/

    conda params.readQC_env

    publishDir "${params.outDir}/db/qc/${params.runID}/", mode: 'copy'

    input:
        path(samplesheet, stageAs: "samplesheet.csv")

    output:
        path('ska_distances.txt'), emit: ska_out
    
    script:
    
    """
    # Make the sample sheet for ska2
        sed 's@,@\t@g' samplesheet.csv \\
            | cut -f2,3,4  \\
            | tail -n +2 \\
            > input.txt

    # Run SKA profiling
        ska build \\
            -o seqs -f input.txt \\
            -k ${params.ska_kmer} \\
            --min-count ${params.ska_min_count} \\
            --proportion-reads ${params.ska_proportion_reads} \\
            --qual-filter ${params.ska_qual_filter} \\
            --min-qual ${params.ska_min_qual} \\
            --threads ${task.cpus}

        ska distance \\
            -o ska_distances.txt \\
            seqs.skf \\
            --min-freq ${params.ska_min_freq} \\
            ${params.ska_distance_args}
    """
}