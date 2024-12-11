process MTBSEQ_SAMPLE_FILTER {
    tag "${runID}"

    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/", mode: 'copy'

    input:
        val(runID)
        val(min_cov)
        val(lineage_pairwise)
        tuple val(sampleID), path(position_table), path(variant_table)
        path(mtbseq_compiled_map_stats)
        path(tbprof_tbdb_res)

    output:
        tuple val(sampleID), path(position_table), path(variant_table), val(lineage), emit: filtered_pairwise_samples
        path "mtbseq.paths.minCov.paths"

    script:
    """
    # Filter samples based on minimum coverage and lineage
    awk -v min_cov="${min_cov}" '\$19 >= min_cov { print \$4 }' ${mtbseq_compiled_map_stats} > samples_with_min_cov.txt

    # Extract lineages for samples with minimum coverage
    grep -f samples_with_min_cov.txt ${tbprof_tbdb_res} | cut -f 1,3 > samples_with_min_cov.lineages.txt

    # Filter samples based on specified lineages
    awk -v lineages="${lineage_pairwise}" '
    BEGIN {split(lineages, lin_arr, ",")}
    {
        for (i in lin_arr) {
            if (\$2 == lin_arr[i]) {
                print \$0
                next
            }
        }
    }
    ' samples_with_min_cov.lineages.txt > filtered_samples.txt

    # Create the output file with filtered samples and their paths
    while IFS=\$'\\t' read -r sample lineage; do
        echo "\${sample}\\t\${position_table}\\t\${variant_table}\\t\${lineage}" >> mtbseq.paths.minCov.paths
    done < filtered_samples.txt

    # Extract lineage for the current sample
    lineage=\$(grep "^${sampleID}\\s" filtered_samples.txt | cut -f2)

    # If lineage is empty, set it to "unknown"
    if [ -z "\$lineage" ]; then
        lineage="unknown"
    fi

    echo "${sampleID}\\t${position_table}\\t${variant_table}\\t\${lineage}" > sample_tuple.txt
    """
}