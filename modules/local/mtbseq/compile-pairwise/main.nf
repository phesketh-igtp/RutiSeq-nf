process MTBSEQ_PAIRWISE_RESULTS {
    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/Clusters/", mode: 'copy'

    input:
        path cluster_files
        path matrix_files

    output:
        path "Master-clusters.csv",     emit: master_clusters
        path "Master-matrixes.csv",     emit: master_matrices

    script:

        def additional_R_script_dir = task.ext.r_script_dir ?: ''

    """

    # Concatenate cluster files
    cat ${cluster_files} > all_clusters.tsv

    Rscript ${additional_R_script_dir}/cluster-results.R \\
                -i all_clusters.tsv \\
                -o Master-cluster.list.tsv

    # Concatenate matrix files
    mkdir corrMatrices/

    for matrix in ${matrix_files}; do

        Rscript ${additional_R_script_dir}/correct-mtbseq-mat.R \\
                -i \${matrix} \\
                -o corrMatrices/\${matrix}.tsv

    done

    Rscript ${additional_R_script_dir}/merge-mtbseq-mat.R \\
                -d corrMatrices/
    
    """

}