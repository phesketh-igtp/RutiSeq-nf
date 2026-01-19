process SNIPPY_SINGLE {

    tag "$sampleID"

    conda params.snippy_env
    
    container { 
        if (workflow.containerEngine == 'singularity') return params.singularity_snp_profiling
        else if (workflow.containerEngine == 'docker') return params.docker_snp_profiling
        else if (workflow.containerEngine == 'apptainer') return params.apptainer_snp_profiling
        else return null
    }
    publishDir "${params.outDir}/db/samples/${sampleID}/snippy/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple val(sampleID), 
            path(fastq_1), 
            path(fastq_2), 
            val(type),
            path(mtbseq_class), 
            path(mtbseq_stats), 
            path(mtbseq_pos), 
            path(mtbseq_vars), 
            path(tbdb_out), 
            path(who_out), 
            path(snippy_out)

    output:
        // tuple for updating the sample ch
        tuple val(sampleID), 
            path(fastq_1), 
            path(fastq_2), 
            val(type),
            path(mtbseq_class), 
            path(mtbseq_stats), 
            path(mtbseq_pos), 
            path(mtbseq_vars),  
            path(tbdb_out),
            path(who_out), 
            path("${sampleID}.vcf"), emit: updated_sample_ch4

        path("${sampleID}.aligned.fa")
        path("${sampleID}.consensus.fa.gz")
        path("${sampleID}.vcf")

    when:
    task.ext.when == null || task.ext.when

    script:
    def additional_args = task.ext.additional_args ?: ''

    """
    # Check if Snippy has already been run for this sample by looking for key output files. Handles situations when the workflow is re-run.
    ## and prevent each single-wf steps from re-running unnecessarily.
        if [[ -f "${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.aligned.fa" \\
                && -f "${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.consensus.fa.gz" \\
                && -f "${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.vcf" \\
                ]]; then

            echo "Snippy results already exist for sample ${sampleID}, skipping Snippy step."
            cp ${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.aligned.fa .
            cp ${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.consensus.fa.gz .
            cp ${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.vcf .

        else
        
        echo -e "Running Snippy for sample ${sampleID}"

            ## determine if we have paired-end or single-end reads

                # Check if we have paired-end reads
                if [[ -f "${fastq_2}" && -s "${fastq_2}" ]]; then
                    # Check if R2 file has actual sequencing data (more than just empty lines)
                    if [[ \$(zcat ${fastq_2} 2>/dev/null | head -n 4 | wc -l) -eq 4 ]]; then
                        echo "Using paired-end mode for ${sampleID}"
                        snippy \\
                            --outdir . \\
                            --force \\
                            --R1 ${fastq_1} \\
                            --R2 ${fastq_2} \\
                            --reference ${params.snippy_reference} \\
                            --mincov ${params.snippy_mincov} \\
                            --minfrac ${params.snippy_minfrac} \\
                            --mapqual ${params.snippy_mapqual} \\
                            --minqual ${params.snippy_minqual} \\
                            --basequal ${params.snippy_basequal} \\
                            --cpus ${task.cpus} \\
                            ${additional_args}
                    else
                        echo "R2 file exists but appears empty, using single-end mode for ${sampleID}"
                        snippy \\
                            --outdir . \\
                            --force \\
                            --se ${fastq_1} \\
                            --reference ${params.snippy_reference} \\
                            --mincov ${params.snippy_mincov} \\
                            --minfrac ${params.snippy_minfrac} \\
                            --mapqual ${params.snippy_mapqual} \\
                            --minqual ${params.snippy_minqual} \\
                            --basequal ${params.snippy_basequal} \\
                            --cpus ${task.cpus} \\
                            ${additional_args}
                    fi
                else
                    echo "Using single-end mode for ${sampleID}"
                    snippy \\
                        --outdir . \\
                        --force \\
                        --se ${fastq_1} \\
                        --reference ${params.snippy_reference} \\
                        --mincov ${params.snippy_mincov} \\
                        --minfrac ${params.snippy_minfrac} \\
                        --mapqual ${params.snippy_mapqual} \\
                        --minqual ${params.snippy_minqual} \\
                        --basequal ${params.snippy_basequal} \\
                        --cpus ${task.cpus} \\
                        ${additional_args}
                fi

            # Rename output files to match expected names
            mv snps.vcf ${sampleID}.vcf
            mv snps.aligned.fa ${sampleID}.aligned.fa
            mv snps.consensus.fa ${sampleID}.consensus.fa

            # Compress consensus file
            gzip --best ${sampleID}.consensus.fa
    
    fi
    """

    stub:
    """
    touch ${sampleID}.aligned.fa
    touch ${sampleID}.consensus.fa
    gzip ${sampleID}.consensus.fa
    touch ${sampleID}.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy: \$(snippy --version 2>&1 | grep -o 'snippy [0-9.]*' | cut -d' ' -f2)
    END_VERSIONS
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