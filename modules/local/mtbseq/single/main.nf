process MTBSEQ_SINGLE {

    tag "$sampleID"

    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/mtbseq").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/mtbseq" : "./modules/local/mtbseq/mtbseq.yml" }

    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data"

    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}", mode: 'mv'

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
    
    """
    # Rename the loaded data to the correct format using the alias/sampleID
    mv ${forward} ${sampleID}_R1.fastq.gz
    mv ${reverse} ${sampleID}_R2.fastq.gz

    echo "Current working directory: \$(pwd)"
    echo "Contents of current directory:"
    ls -l
    echo "Contents of sample file:"
    cat ${sampleID}_sample.txt

    MTBseq --step TBfull --thread ${task.cpus} --window 10

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