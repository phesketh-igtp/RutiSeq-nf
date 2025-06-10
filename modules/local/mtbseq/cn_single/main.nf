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

    tag "$sampleID"

    conda params.mtbseq_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_mtbseq
            else if (workflow.containerEngine == 'docker') return params.docker_mtbseq
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_mtbseq
            else return null
        }
    
    publishDir "${params.outDir}/negative-controls/mtbseq/", mode: 'copy'

    input:
        tuple val(sampleID), 
                path(forward),
                path(reverse)
                
    output:
            path("Classification/${sampleID}.Strain_Classification.tab"), emit: cn_mtbseq_class
            path("Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab"), emit: cn_mtbseq_stats

    script:

        """
        # Rename the reads to the intended naming structure
            if [ "${forward}" != "${sampleID}_R1.fastq.gz" ]; then
                mv ${forward} ${sampleID}_R1.fastq.gz
            fi

            if [ "${reverse}" != "${sampleID}_R2.fastq.gz" ]; then
                mv ${reverse} ${sampleID}_R2.fastq.gz
            fi

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

        mkdir -p Classification/ Statistics/

        # Rename
            if [ -f "Classification/Strain_Classification.tab" ]; then
                mv Classification/Strain_Classification.tab Classification/${sampleID}.Strain_Classification.tab
            else
                echo "" > Classification/${sampleID}.Strain_Classification.tab
            fi

            if [ -f "Statistics/Mapping_and_Variant_Statistics.tab" ]; then
                mv Statistics/Mapping_and_Variant_Statistics.tab Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab
            else 
                echo "" > Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab
            fi

        # Restore read names to original
            if [ "${sampleID}_R1.fastq.gz" != "${forward}" ]; then
                mv ${sampleID}_R1.fastq.gz ${forward}
            fi

            if [ "${sampleID}_R2.fastq.gz" != "${reverse}" ]; then
                mv ${sampleID}_R2.fastq.gz ${reverse}
            fi
        """

}