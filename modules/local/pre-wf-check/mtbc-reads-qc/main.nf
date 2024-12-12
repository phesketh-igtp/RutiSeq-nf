process MTBC_READ_QC {
    
    tag "$sampleID"
    
    conda params.kaiju_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f00cd356ee92f5211e5941beeb4bcab6abfb341e0e5fa7ace8c043406c13381/data'
        } else { 'community.wave.seqera.io/library/kaiju_seqkit:6e4140ab47bd567e' }
    }

    publishDir "${params.outdir}/bbdd/read-qc", mode: 'link'

    input:
        tuple val(sampleID), path(forward), path(reverse)
        path kaiju_names
        path kaiju_nodes
        path kaiju_fmi

    output:
        tuple val(sampleID), 
            path("mtbc_reads/${sampleID}_R1.fastq.gz"), 
            path("mtbc_reads/${sampleID}_R2.fastq.gz"),         emit: mtbc_reads, optional: true
        tuple val(sampleID), path("${sampleID}.qc.out"),          emit: qc_results, optional: true
        path("${sampleID}.kaiju.out"), optional: true
        path("${sampleID}.kaiju_summary.tsv"), optional: true

    script:
    def additional_args_kaiju       = task.ext.additional_args_kaiju ?: ''
    def additional_args_kaiju2table = task.ext.additional_args_kaiju2table ?: ''

    """
    grep 'Mycobacterium tuberculosis ' ${kaiju_names} | cut -f1 | sort | uniq > MTBC.list

    seqkit stats -abT -j ${task.cpus} ${forward} | sed "1d" > tmp.stats.r1.orig
    seqkit stats -abT -j ${task.cpus} ${reverse} | sed "1d" > tmp.stats.r2.orig

    cut -f4 tmp.stats.r1.orig > tmp.r1_num; cut -f17 tmp.stats.r1.orig > tmp.r1_qual
    cut -f4 tmp.stats.r2.orig > tmp.r2_num; cut -f17 tmp.stats.r2.orig > tmp.r2_qual

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

    mkdir -p mtbc_reads/
    seqkit grep -j ${task.cpus} -f tmp.${sampleID}.list ${forward} -o mtbc_reads/${sampleID}_R1.fastq.gz
    seqkit stats -abT -j ${task.cpus} mtbc_reads/${sampleID}_R1.fastq.gz | sed "1d" > tmp.stats.r1.filt

    seqkit grep -j ${task.cpus} -f tmp.${sampleID}.list ${reverse} -o mtbc_reads/${sampleID}_R2.fastq.gz
    seqkit stats -abT -j ${task.cpus} mtbc_reads/${sampleID}_R2.fastq.gz | sed "1d" > tmp.stats.r2.filt

    cut -f4 tmp.stats.r1.filt > tmp.r1_num.filt ; cut -f17 tmp.stats.r1.filt > tmp.r1_qual.filt
    cut -f4 tmp.stats.r2.filt > tmp.r2_num.filt ; cut -f17 tmp.stats.r2.filt > tmp.r2_qual.filt    

    grep 'Mycobacterium' "${sampleID}.kaiju_summary.tsv" | cut -f2 > tmp.MTB.perc
    echo -e ${sampleID} > tmp.sampleID

    # Create the output file of the results
    paste -d '\t' tmp.sampleID \\
                tmp.r1_num \\
                tmp.r1_qual \\
                tmp.r2_num \\
                tmp.r2_qual \\
                tmp.MTB.perc \\
                tmp.r1_num.filt \\
                tmp.r1_qual.filt \\
                tmp.r2_num.filt \\
                tmp.r2_qual.filt \\
                > ${sampleID}.qc.out

    # Clean up temporary files
    rm tmp.*
    """
}