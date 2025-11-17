process SKA_PROFILING {

    conda params.snp_profiling_env

    input:
        path(sylph_classification_report, stageAs: 'sylph_classification.csv')

    output:
        path 'ska_profiling_results.csv'
    
    script:
    
    """
    # Make the sample sheet for ska
        cp ${params.samplesheet} samplesheet.csv
        cut -d ',' -f2,3,4 samplesheet.csv \\
            | sed 's@,@\t@' \\
            > input.txt

    # Run SKA profiling
        ska build \\
            -o seqs -f input.txt \\
            -k ${params.ska_kmer} \\
            --min-count ${params.ska_min_count} \\
            --proportion-reads ${params.ska_proportion_reads} \\\\
            --qual-filter ${params.ska_qual_filter} \\
            --min-qual ${params.ska_min_qual} \\
            --threads ${task.cpus}

        ska distance -o distances.txt \\
            seqs.skf \\
            --min-freq ${params.ska_min_freq} \\
            ${params.ska_distance_args}
    """
}