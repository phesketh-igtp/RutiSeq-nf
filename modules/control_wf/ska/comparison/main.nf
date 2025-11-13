process SKA_MERGE_AND_COMPARISON {
    input:
        path(ska_sketch, stageAs: 'seqs.skf')

    output:
        path 'ska_merge_comparison_results.csv'
    
    shell:
    
    """
        ska distance -o distances.txt .skf
    """
}