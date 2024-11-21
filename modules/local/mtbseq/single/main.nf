process MTBSEQ_SINGLE {

    conda "bioconda::mtbseq=1.1.0"

    container "oras://community.wave.seqera.io/library/mtbseq:1.1.0--d0a0774ca038b4d3"

    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}", mode: 'link'

    input:
    tuple val(sampleID), path(forward, name: "${sampleID}_R1.fastq.gz"), path(reverse, name: "${sampleID}_R2.fastq.gz")

    output:
    path "Bam/",
        path "Bam/${sampleID}.bam",
        path "Bam/${sampleID}.bam.bai",
        path "Bam/${sampleID}.bamlog",
    path "Position_Tables", emit: position_tables_dir
        path "Position_Tables/${sampleID}.gatk_position_table.tab", emit: position_tables
    path "Classification", emit: classification_dir
        path "Classification/Strain_Classification.tab", emit: classification
    path "Statistics", emit: statistics_dir
        path "Statistics/Mapping_and_Variant_Statistics.tab", emit: statistics
    path "Called/", emit: position_variants, emit: called_dir
    path "Called/*gatk_position_variants*.tab", emit: position_variants


    script:
    
    """
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
        --thread ${task.cpus} \
        --window 10
        ${args ?: ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mtbseq: \$(MTBseq --version 2>&1 | sed 's/^.*MTBseq //; s/ .*\$//')
    END_VERSIONS
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