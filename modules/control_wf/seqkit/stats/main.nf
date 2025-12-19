process READS_STATS {

/*
    @author: Poppy J hesketh Best
    @date: 2025-11-17
    @version: v1.1.0
    @description:
        Produces read statistics of the input fasta files with seqkit.
    @changelog:
        v1.0.0-2025-11-17: Functioning module created.
        v1.1.0-2025-12-19: Modified to handle when the samlesheet has an empty line
*/

    conda params.readQC_env

    publishDir "${params.outDir}/db/qc/${params.runID}/"

    input:
        path(samplesheet, stageAs: 'samplesheet.csv')

    output:
    path("seqkit_stats.txt"), emit: reads_stats_out

    script:
    """
    mkdir -p reads/
    
    while IFS=',' read -r identities samples fastq_1 fastq_2 type; do
        if [[ -n \$fastq_1 ]]; then
            ln -s \$fastq_1 reads/
        fi

        if [[ -n \$fastq_2 ]]; then
            ln -s \$fastq_2 reads/
        fi
    done < <(tail -n +2 samplesheet.csv)

    # Generate fastq stats with seqkit
    seqkit stat -bT reads/*.fastq.gz \\
        > seqkit_stats.txt
    """
}