process READS_STATS {

/*
    @author: Poppy J hesketh Best
    @date: 2025-11-17
    @version: v1.0.0
    @description:
        Produces read statistics of the input fasta files with seqkit.
    @changelog:
        v1.0.0-2025-11-17: Functioning module created.
*/

    conda params.readQC_env

    storeDir "${params.outDir}/db/qc/${params.runID}/"

    input:
        path(samplesheet, stageAs: 'samplesheet.csv')

    output:
    path("seqkit_stats.txt"), emit: reads_stats_out

    script:
    """
    mkdir -p reads/
    
    # Create symbolic links to the fastq files
    while IFS=',' read -r identities samples fastq_1 fastq_2 type; do
        ln -s \${fastq_1} reads/
        ln -s \${fastq_2} reads/
    done < <(tail -n +2 samplesheet.csv)

    # Generate fastq stats with seqkit
    seqkit stat -bT reads/*.fastq.gz \\
        > seqkit_stats.txt
    """
}