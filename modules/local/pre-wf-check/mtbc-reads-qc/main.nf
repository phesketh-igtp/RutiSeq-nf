process MTBC_READ_QC {
    tag "$sampleID"
    
    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/kaiju").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/kaiju" : "../modules/local/pre-wf-check/mtbc-reads-qc/kaiju.yml" }
    
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/09/09c601d67c915ec87dbcbf5b8a12e5cc15d6bc36523377b5f3103c3ec2dc0c50/data'
    
    cpus { Math.min(16 * task.attempt, 24) }
        
    input:
        tuple val(sampleID), path(forward), path(reverse)

    output:
        tuple val(sampleID), path("${sampleID}_R1.fastq.gz"), path("${sampleID}_R2.fastq.gz"), emit: mtbc_reads
        path("${sampleID}.qc.out"), emit: qc_out

    script:
    """
    grep 'Mycobacterium tuberculosis ' ${params.kaiju_names} | cut -f1 | sort | uniq > MTBC.list
    echo ${sampleID} > sampleID

    seqkit stats -abT -j ${task.cpus} ${forward} | sed '1d' > stats.r1
    seqkit stats -abT -j ${task.cpus} ${reverse} | sed '1d' > stats.r2
    
    cut -f4 stats.r1 > r1_num; cut -f17 stats.r1 > r1_qual
    cut -f4 stats.r2 > r2_num; cut -f17 stats.r2 > r2_qual

    kaiju -t ${params.kaiju_nodes} -f ${params.kaiju_fmi} \
        -i ${forward} -j ${reverse} \
        -z ${task.cpus} -o ${sampleID}.kaiju.out
    
    kaiju2table -t "${params.kaiju_nodes}" -n "${params.kaiju_names}" \
                -r genus -m 1.0 -o "${sampleID}.kaiju_summary.tsv" \
                "${sampleID}.kaiju.out"
    
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