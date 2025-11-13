process READS_STATS {

    conda "${params.projectDir}/env/conda/snp_profiling_env.yml"

    input:
        path(samplesheet, stageAs: 'samplesheet.csv')

    output:
    path("seqkit_stats.txt"), emit: reads_stats_out

    script:
    """
    # Create symbolic links to the fastq files
    while read -r identities samples fastq_1 fastq_2 type; do
        ln -s \${fastq_1} .
        ln -s \${fastq_2} .
    done < <(tail -n +2 samplesheet.csv)

    # Generate fastq stats with seqkit
    seqkit stat -bT *.fasta.gz \\
        > seqkit_stats.txt
    """
}