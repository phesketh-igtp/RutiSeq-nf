process MTBSEQ_SINGLE {

    conda "bioconda::mtbseq=1.1.0"

    container "oras://community.wave.seqera.io/library/mtbseq:1.1.0--d0a0774ca038b4d3"

    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}", mode: 'link'

    input:
    tuple val(sampleID), path(forward), path(reverse)

    output:
    path "Amend", emit: called
    path "Position_Tables", emit: position_tables_dir
    path "Classification", emit: classification_dir
    path "Statistics", emit: statistics_dir
    path "Statistics/Mapping_and_Variant_Statistics.tab", emit: statistics
    path "Classification/Strain_Classification.tab", emit: classification
    path "Called/*gatk_position_variants*.tab", emit: position_variants
    path "Position_Tables/*.gatk_position_table.tab", emit: position_tables

    script:
    
    """

    # Create symbolic links for input files using full paths
    ln -s ${forward.toRealPath()} ${sampleID}_R1.fastq.gz
    ln -s ${reverse.toRealPath()} ${sampleID}_R2.fastq.gz

    unlink ${forward.toRealPath()}
    unlink ${reverse.toRealPath()}

    # Verify that the symlinks were created successfully
    ls -l ${sampleID}_R1.fastq.gz ${sampleID}_R2.fastq.gz

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