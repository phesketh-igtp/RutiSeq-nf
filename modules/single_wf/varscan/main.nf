process VARSCAN_SINGLE {

    tag "$sampleID"

    conda params.varscan_env
    
    container { 
        if (workflow.containerEngine == 'singularity') return params.singularity_snp_profiling
        else if (workflow.containerEngine == 'docker') return params.docker_snp_profiling
        else if (workflow.containerEngine == 'apptainer') return params.apptainer_snp_profiling
        else return null
    }
    publishDir "${params.outDir}/db/samples/${sampleID}/varscan/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple val(sampleID), 
            path(fastq_1), 
            path(fastq_2), 
            path(bam)

    output:
        // tuple for updating the sample ch
        path("${sampleID}.varscan.vcf")
        path("${sampleID}.varscan.consensus.fa")

    when:
    task.ext.when == null || task.ext.when

    script:
    def additional_args = task.ext.additional_args ?: ''

    """
    # Check if Snippy has already been run for this sample by looking for key output files. Handles situations when the workflow is re-run.
    ## and prevent each single-wf steps from re-running unnecessarily.
    if [[ -f "${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.bam" \\
        && -f "${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.varscan.vcf" \\
        && -f "${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.varscan.consensus.fa"
    ]]; then

        echo "Snippy results already exist for sample ${sampleID}, skipping Snippy step."
        ln -s ${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.bam .
        ln -s ${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.varscan.vcf
        ln -s ${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.varscan.consensus.fa .

    else
        
        echo -e "Running VarScan2 for sample ${sampleID}"
        samtools mpileup \\
            -q ${params.varscan_min_avg_qual} \\
            -f ${params.snippy_reference_fa} \\
            ${sampleID}.bam \\
            -o ${sampleID}.mpileup

        varscan mpileup2cns \\
            ${sampleID}.mpileup \\
            --min-avg-qual ${params.varscan_min_avg_qual} \\
            --min-coverage ${params.varscan_min_coverage} \\
            --variants \\
            --output-vcf 1 \\
            --strand-filter 0 \\
            > ${sampleID}.varscan.vcf
        sed -i "s@Sample1@${sampleID}@g" ${sampleID}.varscan.vcf

    # Compress and index
        bcftools convert \\
            -Oz -o ${sampleID}.varscan.vcf.gz \\
            ${sampleID}.varscan.vcf
        bcftools index -f ${sampleID}.varscan.vcf.gz

        # Call the consensus:
        bcftools consensus \\
            --sample '${sampleID}' \\
            --fasta-ref ${params.snippy_reference_fa} \\
            --output ${sampleID}.varscan.consensus.fa \\
            ${sampleID}.varscan.vcf.gz
            
        # Compress consensus file
        gzip --best ${sampleID}.varscan.consensus.fa

    fi
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2026-01-19
@version: 1.0
@description: 
    This process runs the SNP profiling step using a reference GBK file.
    It produces files that can used for snippy-core processing
    To generate large phylogenetic trees. Produces a consensus sequence.
@changelog:
    v1.0.0-2025-04-01: Initial version, using MTBSeq mpileup output
    v2.0.0-2025-06-10: Change from using MTBSeq outputs to using Snippy for SNP profiling
                        Added support for single-end reads
                        Extended available parameters for Snippy
    v2.1.0-2026-01-19: Updated to check for existing results to avoid re-running
*/