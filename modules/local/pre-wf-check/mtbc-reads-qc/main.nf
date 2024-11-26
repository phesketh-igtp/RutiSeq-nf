process MTBC_READ_QC {

    tag "$sampleID"
    
    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/kaiju").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/kaiju" : "../modules/local/pre-wf-check/mtbc-reads-qc/kaiju.yml" }
    
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/09/09c601d67c915ec87dbcbf5b8a12e5cc15d6bc36523377b5f3103c3ec2dc0c50/data'
        
    input:
        tuple val(sampleID), path(forward), path(reverse), val(db_done)

    output:
        tuple val(sampleID), path("mtbc-reads/${sampleID}_R1.fastq.gz"), path("mtbc-reads/${sampleID}_R2.fastq.gz"), emit: mtbc_reads
        path("samples/${sampleID}.qc.out"), emit: qc_out

    script:
    """
    grep 'Mycobacterium tuberculosis ' "${params.kaiju_names}" | cut -f1 | sort | uniq > .tmp.MTBC.list

    seqkit stats -abT -j "${params.cpu}" "${forward}" | sed '1d' > .${sampleID}.tmp-r1
    seqkit stats -abT -j "${params.cpu}" "${reverse}" | sed '1d' > .${sampleID}.tmp-r2
    R1_num=\$(cut -f4 .${sampleID}.tmp-r1)
    R1_qual=\$(cut -f17 .${sampleID}.tmp-r1)
    R2_num=\$(cut -f4 .${sampleID}.tmp-r2)
    R2_qual=\$(cut -f17 .${sampleID}.tmp-r2)

    kaiju -t "${params.kaiju_nodes}" -f "${params.kaiju_fmi}" \
            -i "${forward}" -j "${reverse}" \
            -z "${params.cpu}" -o "${sampleID}.kaiju.out"
    
    kaiju2table -t "${params.kaiju_nodes}" -n "${params.kaiju_names}" \
                -r genus -m 1.0 -o "${sampleID}.kaiju_summary.tsv" \
                "${sampleID}.kaiju.out"
    
    #MTB_PERC=\$(grep 'Mycobacterium' "${sampleID}.kaiju_summary.tsv" | sed '1d' | head -1 | cut -f2)

    mkdir -p mtbc-reads/
    #echo -e "\${sampleID}\t\${R1_num}\t\${R1_qual}\t\${R2_num}\t\${R2_qual}\t\${MTB_PERC}" > samples/${sampleID}.qc.out

    grep -f .tmp.MTBC.list "${sampleID}.kaiju.out" | cut -f2 > .tmp.${sampleID}.list
    seqkit grep -j "${params.cpu}" -f .tmp.${sampleID}.list "${forward}" -o mtbc-reads/${sampleID}_R1.fastq.gz
    seqkit grep -j "${params.cpu}" -f .tmp.${sampleID}.list "${reverse}" -o mtbc-reads/${sampleID}_R2.fastq.gz

    # Clean up temporary files
    rm .tmp.* .${sampleID}.tmp-r1 .${sampleID}.tmp-r2
    """
}