process FETCH_SRA {

    tag "$sampleID"
    
    // conda directive is ignored when using containers
    conda 'bioconda::sra-tools=3.0.10'
    
    input:
        tuple val(sampleID), 
            path(forward), path(reverse), 
            val(type), path(mtbseq_class),
            path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
            path(tbdb_out), path(who_out), path(mtbseq_vcf)
    
    output:
    tuple val(sampleID), 
        path("${sampleID}_R1.fastq.gz"), 
        path("${sampleID}_R2.fastq.gz", optional: true),
        val("sample"), 
        path(mtbseq_class),
        path(mtbseq_stats),
        path(mtbseq_pos),
        path(mtbseq_vars), 
        path(tbdb_out),
        path(who_out),
        path(mtbseq_vcf), emit: fetched_sra_tuple
        
    script:
    """
    # Download FASTQ files
    fastq-dump --split-files --gzip --outdir . ${forward}
    
    # Handle single-end vs paired-end data
    if [[ -f "${forward}_2.fastq.gz" ]]; then
        echo "Paired-end data detected for ${forward}"
        mv ${forward}_1.fastq.gz ${forward}_R1.fastq.gz
        mv ${sampleID}_2.fastq.gz ${sampleID}_R2.fastq.gz
    else
        echo "Single-end data detected for ${forward}"
    fi
    """
}


/*
@author: Poppy J Hesketh Best
@date: 2025-11-13
@version: v1.0.0
@description:
    The download of SRA reads using fastq-dump, and creation of a empty
    fastq_2  if the data is SE reads, to satisfy the output format.
@changelog:
    2025-11-13.v1.0.0: Initial Version
*/
