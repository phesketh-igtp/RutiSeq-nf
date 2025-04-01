process CN_MTBSEQ_SINGLE {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 0.1
    @description: 
        This module runs MTBseq on a single sample. It is designed to be used in the context 
        of the negative control workflow. It takes a tuple of sampleID, forward read file, 
        and reverse read file as input. The output is the MTBseq classification and 
        statistics files.
*/

    tag "${sampleID}"

    conda params.mtbseq_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_mtbseq
            else if (workflow.containerEngine == 'docker') return params.docker_mtbseq
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_mtbseq
            else return null
        }
    
    publishDir "${params.outdir}/bbdd/negative-controls/mtbseq/${sampleID}", mode: 'copy'

    input:
        tuple val(sampleID), 
                path(forward),
                path(reverse)
                
    output:
        path("Classification/${sampleID}.Strain_Classification.tab")
        path("Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab")

    script:

        """
        # Rename the reads to the intended naming structure
            mv ${forward} ${sampleID}_R1.fastq.gz
            mv ${reverse} ${sampleID}_R2.fastq.gz

        # Run MTBseq for a single sample
            MTBseq --step TBfull \\
                    --thread ${task.cpus} \\
                    --minbqual ${params.mtbseq_minbqual} \\
                    --mincovf ${params.mtbseq_mincovf} \\
                    --mincovr ${params.mtbseq_mincovr} \\
                    --minphred20 ${params.mtbseq_minphred20} \\
                    --minfreq ${params.mtbseq_minfreq} \\
                    --unambig ${params.mtbseq_unambig} \\
                    --window ${params.mtbseq_window} \\
                    1>>.command.out \\
                    2>>.command.err || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

        # Renamecp Classification/Strain_Classification.tab Classification/${sampleID}.Strain_Classification.tab
            if [ -f "Classification/Strain_Classification.tab" ]; then
                mv Classification/Strain_Classification.tab Classification/${sampleID}.Strain_Classification.tab
            fi

            if [ -f "Statistics/Mapping_and_Variant_Statistics.tab" ]; then
                mv Statistics/Mapping_and_Variant_Statistics.tab Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab
            fi

        # restore the symbolic link names
            mv ${sampleID}_R1.fastq.gz ${forward}; mv ${sampleID}_R2.fastq.gz ${reverse}

        # Always exit with status 0 to prevent pipeline failure
            exit 0
        """

}