process CN_READ_TAXONOMY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-04
    @version: 0.2
    @description: 
        Run KAIJU on the reads and get read taxonomy and statistics of the reads. Outputs of this modules 
        are intended to be combines into a single file for each sample, and then concatenated into a
        single file for all samples for a particular run.

    @updates: v0.2 - created new paths for the results `Classification/ Statistics/`
*/
    
    tag "${sampleID}"

    conda params.taxonomy_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f00cd356ee92f5211e5941beeb4bcab6abfb341e0e5fa7ace8c043406c13381/data'
        } else { 'community.wave.seqera.io/library/kaiju_seqkit:6e4140ab47bd567e' }
    }

    publishDir "${params.outdir}/bbdd/negative-controls/", mode: 'link'

    input:
        tuple val(sampleID), 
                path(forward), 
                path(reverse)

    output:
        path("Classification/${sampleID}.k2.report"),   emit: cn_k2_report
        path("Statistics/${sampleID}.stats.tsv"),       emit: cn_stats
        path("Classification/${sampleID}.k2.output.gz")

    script:

        """
        mkdir -p Classification/ Statistics/

        # Change read names to be used in kraken2
            mv ${forward} ${sampleID}_R1.fastq.gz
            mv ${reverse} ${sampleID}_R2.fastq.gz

        # Run kraken2
            kraken2 \\
                --threads ${task.cpus} \\
                --db ${params.kraken_db_path} \\
                --memory-mapping \\
                --report Classification/${sampleID}.k2.report \\
                --paired ${sampleID}_R1.fastq.gz ${sampleID}_R2.fastq.gz \\
                > Classification/${sampleID}.k2.output

            gzip --best Classification/${sampleID}.k2.output

        # Genreate stats
            seqkit stats -bTa ${sampleID}_R1.fastq.gz | sed 's@.fastq.gz@@g' > tmp.for.tsv
            seqkit stats -bTa ${sampleID}_R2.fastq.gz | sed 's@.fastq.gz@@g' | sed '1d' > tmp.rev.tsv

            cat tmp.for.tsv tmp.rev.tsv > Statistics/${sampleID}.stats.tsv
            rm tmp.for.tsv tmp.rev.tsv
        
        # Add a new column to the results containing the sampleID for later concatenation
            sed -i "s@^@${sampleID}\t@" Classification/${sampleID}.k2.report
            sed -i "s@^@${sampleID}\t@" Statistics/${sampleID}.stats.tsv

        # Restore read names to original
        mv ${sampleID}_R1.fastq.gz ${forward} 
        mv ${sampleID}_R2.fastq.gz ${reverse} 
        """
}