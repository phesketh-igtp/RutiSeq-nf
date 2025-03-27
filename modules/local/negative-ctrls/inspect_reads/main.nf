process CN_READ_TAXONOMY {
    
    tag "$sampleID"

    conda params.taxonomy_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f00cd356ee92f5211e5941beeb4bcab6abfb341e0e5fa7ace8c043406c13381/data'
        } else { 'community.wave.seqera.io/library/kaiju_seqkit:6e4140ab47bd567e' }
    }

    publishDir "${params.outdir}/bbdd/negative-controls/results/", mode: 'link'

    input:
        tuple val(sampleID), 
            path(forward), 
            path(reverse),
            path(qc_results)

    output:
        path("${sampleID}.k2.report"),     emit: cn_k2_report
        path("${sampleID}.k2.raw.tsv")
        path("${sampleID}.stats.tsv"),  emit: cn_stats

    script:

        """
        # Change read names to be used in kraken2
        mv ${forward} ${sampleID}_R1.fastq.gz
        mv ${reverse} ${sampleID}_R2.fastq.gz

        kraken2 \\
            --threads ${task.cpus} \\
            --db ${params.kraken_db_path} \\
            --use-names \\
            --use-mpa-style \\
            --memory-mapping \\
            --report ${sampleID}.k2.report \\
            --paired ${sampleID}_R1.fastq.gz ${sampleID}_R2.fastq.gz \\
            > ${sampleID}.k2.raw.tsv

        seqkit stats -bTa ${sampleID}_R1.fastq.gz > tmp.for.tsv
        seqkit stats -bTa ${sampleID}_R2.fastq.gz | sed '1d' > tmp.rev.tsv

        cat tmp.for.tsv tmp.rev.tsv > ${sampleID}.stats.tsv

        # Restore read names to original
        mv ${sampleID}_R1.fastq.gz ${forward} 
        mv ${sampleID}_R2.fastq.gz ${reverse} 
        """
}