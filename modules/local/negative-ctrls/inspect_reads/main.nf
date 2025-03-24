process CN_READ_TAXONOMY {
    
    tag "$sampleID"

    conda params.tax_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f00cd356ee92f5211e5941beeb4bcab6abfb341e0e5fa7ace8c043406c13381/data'
        } else { 'community.wave.seqera.io/library/kaiju_seqkit:6e4140ab47bd567e' }
    }

    publishDir "${params.outdir}/bbdd/negative-controls/results/", mode: 'link'

    input:
        tuple val(sampleID), path(forward), path(reverse)

    output:
        path("${sampleID}.qc.out"),    emit: cn_qc_results
        path("${sampleID}.kaiju.out")
        path("${sampleID}.kaiju_summary.tsv"), emit: cn_kaiju_results

    script:

        """
        kraken2 --paired ${forward} ${reverse} \\
            --output ${sampleID}.k2.out \\
            --log ${sampleID}.k2.log \\
            --minimum-base-quality 30 \\
            --use-names \\
            --use-mpa-style \\
            --memory-mapping \\
            --threads ${task.cpus} \\
            --db ${params.kraken_db_path}
        """
}