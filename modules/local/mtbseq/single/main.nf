process MTBSEQ_SINGLE {
    tag "$sampleID"

    conda "bioconda::mtbseq=1.1.0"

    container "oras://community.wave.seqera.io/library/mtbseq:1.1.0--d0a0774ca038b4d3"

    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}", mode: 'link'

    input:
    tuple val(sampleID), path(forward), path(reverse)

    output:
    path "Bam/", emit: bam_dir
    path "Bam/${sampleID}.bam", emit: bam
    path "Bam/${sampleID}.bam.bai", emit: bam_index
    path "Bam/${sampleID}.bamlog", emit: bamlog
    path "Position_Tables", emit: position_tables_dir
    path "Position_Tables/${sampleID}.gatk_position_table.tab", emit: position_tables
    path "Classification", emit: classification_dir
    path "Classification/Strain_Classification.tab", emit: classification
    path "Statistics", emit: statistics_dir
    path "Statistics/Mapping_and_Variant_Statistics.tab", emit: statistics
    path "Called/", emit: called_dir
    path "Called/*gatk_position_variants*.tab", emit: position_variants

    script:
    def threads = params.mtbseq_threads ?: task.cpus
    
    """
    # Rename the loaded data to the correct format using the alias/sampleID
    mv ${forward} ${sampleID}_R1.fastq.gz
    mv ${reverse} ${sampleID}_R2.fastq.gz

    # Create sample file
    echo "${sampleID} ${sampleID}_R1.fastq.gz ${sampleID}_R2.fastq.gz" > ${sampleID}_sample.txt

    echo "Current working directory: \$(pwd)"
    echo "Contents of current directory:"
    ls -l
    echo "Contents of sample file:"
    cat ${sampleID}_sample.txt

    MTBseq --step TBfull \
        --samples ${sampleID}_sample.txt \
        --project ${sampleID} \
        --thread ${threads} \
        --window 10

    """

    stub:
    """
    mkdir -p Amend Position_Tables Classification Statistics Called
    touch Statistics/Mapping_and_Variant_Statistics.tab
    touch Classification/Strain_Classification.tab
    touch Called/gatk_position_variants.tab
    touch Position_Tables/gatk_position_table.tab
    """
}