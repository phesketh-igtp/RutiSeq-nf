process CN_MTBSEQ_SINGLE {

    tag "$sampleID"

    conda params.mtbseq_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_mtbseq
            else if (workflow.containerEngine == 'docker') return params.docker_mtbseq
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_mtbseq
            else return null
        }
    
    publishDir "${params.outdir}/bbdd/mtbseq/negative-controls/samples/${sampleID}", mode: 'copy'

    input:
        tuple val(sampleID), 
                    path(forward, stageAs: "${sampleID}_R1.fastq.gz"),
                    path(reverse, stageAs: "${sampleID}_R2.fastq.gz")
                
    output:
        path("Called/*")
        path("Classification/${sampleID}.Strain_Classification.tab")
        path("Mpileup/*")
        path("Position_Tables/*")
        path("Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab")
        path("Bam/*")
        path("GATK_Bam/*")

    script:

        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        # Run MTBseq for a single sample
            MTBseq --step TBfull \\
                --thread        ${task.cpus} \\
                --minbqual      ${params.mtbseq_minbqual} \\
                --mincovf       ${params.mtbseq_mincovf} \\
                --mincovr       ${params.mtbseq_mincovr} \\
                --minphred20    ${params.mtbseq_minphred20} \\
                --minfreq       ${params.mtbseq_minfreq} \\
                --unambig       ${params.mtbseq_unambig} \\
                --window        ${params.mtbseq_window} \\
                ${additional_args} \\
                1>>.command.out \\
                2>>.command.err || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

        # Rename the stats and class outputs to have unique names
            mv Classification/Strain_Classification.tab Classification/${sampleID}.Strain_Classification.tab
            mv Statistics/Mapping_and_Variant_Statistics.tab Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab
        """

}