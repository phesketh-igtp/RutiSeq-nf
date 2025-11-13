process FETCH_SRA {

    tag "$accession"
    
    // conda directive is ignored when using containers
    conda 'bioconda::sra-tools=3.0.10'
    
    input:
        tuple val(accession), 
            path(forward), path(reverse), 
            val(type),  path(mtbseq_class),
            path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
            path(tbdb_out), path(who_out), path(mtbseq_vcf)
    
    output:
    tuple val(accession), 
        path("${accession}_1.fastq.gz"), 
        path("${accession}_2.fastq.gz"),
        val(type),  path(mtbseq_class),
        path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
        path(tbdb_out), path(who_out), path(mtbseq_vcf), emit: fetch_fastq_tuple
        
    script:
    """
    # Download FASTQ files
    fastq-dump --split-files --gzip --outdir . ${accession}
    
    # Handle single-end vs paired-end data
    if [[ -f "${accession}_2.fastq.gz" ]]; then
        echo "Paired-end data detected for ${accession}"
    else
        echo "Single-end data detected for ${accession}, creating empty R2"
        # Create empty file for consistency
        echo -n | gzip > ${accession}_2.fastq.gz
    fi
    """
}