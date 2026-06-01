process MTBSEQ_SINGLE {

    tag "$sampleID"

    conda params.mtbseq_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_mtbseq
            else if (workflow.containerEngine == 'docker') return params.docker_mtbseq
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_mtbseq
            else return null
        }
    
    publishDir "${params.outDir}/db/samples/${sampleID}/mtbseq/", 
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
                path(mtbseq_vcf)
                
    output:
        path("Called/*")
        path("Classification/${sampleID}.Strain_Classification.tab")
        path("Position_Tables/*")
        path("Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab")

        // tuple for updating the sample_ch
        tuple val(sampleID), 
                path(fastq_1), 
                path(fastq_2),
                val(type),
                path("Classification/${sampleID}.Strain_Classification.tab"), 
                path("Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab"), 
                path("Position_Tables/${sampleID}.gatk_position_table.tab"), 
                path("Called/${sampleID}.gatk_position_variants_*.tab"),  
                path(tbdb_out), 
                path(who_out), 
                path(mtbseq_vcf), emit: updated_sample_ch3

    script:

        """
        # Check if MTBSeq has already been run for this sample by looking for key output files. Handles situations when the workflow is re-run.
        ## and prevent each single-wf steps from re-running unnecessarily.
        if [[ -f "${params.outDir}/db/samples/${sampleID}/mtbseq/Classification/${sampleID}.Strain_Classification.tab" \\
            && -f "${params.outDir}/db/samples/${sampleID}/mtbseq/Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab" \\
            && -f "${params.outDir}/db/samples/${sampleID}/mtbseq/Position_Tables/${sampleID}.gatk_position_table.tab" \\
                ]]; then

            echo "MTBSeq results already exist for sample ${sampleID}, skipping MTBSeq step..."
            mkdir -p Called Position_Tables Classification Statistics
            ln -s ${params.outDir}/db/samples/${sampleID}/mtbseq/Classification/${sampleID}.Strain_Classification.tab Classification/
            ln -s  ${params.outDir}/db/samples/${sampleID}/mtbseq/Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab Statistics/
            ln -s  ${params.outDir}/db/samples/${sampleID}/mtbseq/Position_Tables/${sampleID}.gatk_position_table.tab Position_Tables/
            ln -s  ${params.outDir}/db/samples/${sampleID}/mtbseq/Called/${sampleID}.gatk_position_variants_*.tab Called/

        else

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
                ${params.mtbseq_args} \\
                1>>.command.out \\
                2>>.command.err || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

        # Rename the stats and class outputs to have unique names
        ## this prevent clashes later on
            cat Classification/Strain_Classification.tab > Classification/${sampleID}.Strain_Classification.tab
            cat Statistics/Mapping_and_Variant_Statistics.tab > Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab

        

        fi \\
        1>>.command.out \\
        2>>.command.err || true
        """

}

/*
@author: Poppy J Hesketh Best
@date: 2026-01-19
@version: 2.1.0
@description: 
    In this module MTBSeq is run for a single genome up to the strain classifications step,
    and will then stop as there are no other genomes to compare with. To start with the forward/reverse
    reads from the MTBC fitlered have a different naming convention (${sampleID}_mtbc_R1.fastq.gz), so 
    the first step is return the naming convention back to the expected name for the MTBeq outputs - 
    I tried stageAs: as some point to circumvent this issue, but it caused other problems since
    the files needs to have the ${sampleID} in the name for the next steps, so I reverted back to
    just moving the files into a name structure. This is a bit of a hack, but it works for now.
        TODO: revisit this issue.
@chagelog
    v1.0.0-2025-04-01: Initial version
    v1.0.1-2025-04-04: Added more comments and description
    v2.0.0-2025-06-10: Added support for different MTBSeq referencess
    v2.1.0-2026-01-19: Updated to check for existing results to avoid re-running
*/