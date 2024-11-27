process MTBC_READ_QC {
    tag "$sampleID"
    
    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/kaiju").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/kaiju" : "../modules/local/pre-wf-check/mtbc-reads-qc/kaiju.yml" }
    
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f00cd356ee92f5211e5941beeb4bcab6abfb341e0e5fa7ace8c043406c13381/data'

    publishDir "${params.outdir}/bbdd/read-qc", mode: 'link'

    input:
        tuple val(sampleID), path(forward), path(reverse)
        path kaiju_names
        path kaiju_nodes
        path kaiju_fmi

    output:
        tuple val(sampleID), path("${sampleID}_R1.fastq.gz"), path("${sampleID}_R2.fastq.gz"), emit: mtbc_reads
        tuple val(sampleID), path("${sampleID}_S*_L001_R1_001.fastq.gz"), path("${sampleID}_S*_L001_R2_001.fastq.gz"), emit: original_reads
        path("${sampleID}.qc.out"), emit: qc_out

    script:

    def additional_args_kaiju = task.ext.additional_args_kaiju ?: '' // defined in the nextflow.config file
    def additional_args_kaiju2table = task.ext.additional_args_kaiju2table ?: '' // defined in the nextflow.config file

    """
    grep 'Mycobacterium tuberculosis ' ${kaiju_names} | cut -f1 | sort | uniq > MTBC.list
    echo ${sampleID} > sampleID

    seqkit stats -abT -j ${task.cpus} ${forward} | sed '1d' > stats.r1
    seqkit stats -abT -j ${task.cpus} ${reverse} | sed '1d' > stats.r2

    cut -f4 stats.r1 > r1_num; cut -f17 stats.r1 > r1_qual
    cut -f4 stats.r2 > r2_num; cut -f17 stats.r2 > r2_qual

    kaiju -t ${kaiju_nodes} -f ${kaiju_fmi} \\
            -i ${forward} -j ${reverse} \\
            -z ${task.cpus} \\
            ${additional_args_kaiju} \\
            -o ${sampleID}.kaiju.out

    kaiju2table -t "${kaiju_nodes}" -n ${kaiju_names} \\
                ${additional_args_kaiju2table} \\
                -o ${sampleID}.kaiju_summary.tsv \\
                ${sampleID}.kaiju.out

    grep -f MTBC.list "${sampleID}.kaiju.out" | cut -f2 > tmp.${sampleID}.list

    seqkit grep -j ${task.cpus} -f tmp.${sampleID}.list ${forward} -o ${sampleID}_R1.fastq.gz
    seqkit grep -j ${task.cpus} -f tmp.${sampleID}.list ${reverse} -o ${sampleID}_R2.fastq.gz

    # Create the output file of the results
    grep 'Mycobacterium' "${sampleID}.kaiju_summary.tsv" | sed '1d' | head -1 | cut -f2 > MTB.perc
    paste -d ',' sampleID r1_num r1_qual r2_num r2_qual MTB.perc > ${sampleID}.qc.out

    # Clean up temporary files
    rm tmp.* stats.r1 stats.r2 MTBC.list r1_num r1_qual r2_num r2_qual sampleID MTB.perc
    
    """


        
}