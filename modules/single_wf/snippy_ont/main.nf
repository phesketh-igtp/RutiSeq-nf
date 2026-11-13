process SNIPPY_SINGLE_ONT {

    tag "$sampleID"

    conda params.snippy_ng_env
    
    container { 
        if (workflow.containerEngine == 'singularity') return params.singularity_snp_profiling
        else if (workflow.containerEngine == 'docker') return params.docker_snp_profiling
        else if (workflow.containerEngine == 'apptainer') return params.apptainer_snp_profiling
        else return null
    }
    publishDir "${params.outDir}/db/samples/${sampleID}/snippy-ng/", 
        mode: 'copy',
        overwrite: true,
        pattern: '{*.tab,*.txt,*.tsv,*.vcf.gz,*.pass.vcf}'

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
            path(tbdb_out),
            path(who_out), 
            path("${sampleID}.pass.vcf"), emit: updated_sample_ch4

        path("CITATIONS.txt")
        path("FILES.tsv")
        path("${sampleID}.all.vcf.gz")
        path("${sampleID}.all.vcf.summary.tsv")
        path("${sampleID}.pass.vcf")
        path("${sampleID}.qc.tsv")
        path("${sampleID}.tab")
        path("LOG.txt")
        path("PIPELINE.txt")

    when:
    task.ext.when == null || task.ext.when

    script:
    def additional_args = task.ext.additional_args ?: ''

    """
    # Check if Snippy has already been run for this sample by looking for key output files. Handles situations when the workflow is re-run.
    ## and prevent each single-wf steps from re-running unnecessarily.
    if [[ -f "${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.all.vcf.summary.tsv" \\
        && -f "${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.all.vcf.gz" \\
        && -f "${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.pass.vcf" \\
        && -f "${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.qc.tsv" \\
        && -f "${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.tab" \\

        ]]; then
            echo "Snippy results already exist for sample ${sampleID}, skipping Snippy step."
            ln -s ${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.all.vcf.summary.tsv .
            ln -s ${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.all.vcf.gz .
            ln -s ${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.pass.vcf .
            ln -s ${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.qc.tsv .
            ln -s ${params.outDir}/db/samples/${sampleID}/snippy-ng/${sampleID}.tab .

    else
        
        echo -e "Running Snippy for sample ${sampleID}"
        
        echo "Using single-end mode for ${sampleID}"
            snippy long \\
                --outdir . \\
                --force \\
                --se ${fastq_1} \\
                --reference ${params.snippy_reference_gbk} \\
                --cpus ${task.cpus} \\
                --prefix ${sampleID} \\
                ${additional_args}

        

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