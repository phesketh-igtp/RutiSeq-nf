process MTBSEQ_SINGLE {

    /*
        In this module MTBSeq is run for a single genome up till the strain classifications step,
            and will then stop as there are no other genomes to compare with. To start with the forward/reverse
            reads from the MTBC fitlered have a different naming convention (${sampleID}_mtbc_R1.fastq.gz), so 
            the first step is return the naming convention back to the expected name for the MTBeq outputs - 
            I am certain that I tried stageAs: as some point to circumvent this issue, but it caused other problems.
            TODO: revisit this issue.
        Like the previous modules, the input for this module is the paths to all the files
            needed for a sample to proceed into the PAIRWISE_WF(), which litters the publish directory
            these excess files and are removed at the very end. This was originally done with a IF argument
            to only remove them if the files did not exist - this prevents nextflow complaining of they 
            were not generated - but this created some sybtax errors, so I changed it to touching the 
            files then removing them manually. TODO: might be to revisit this at a later date.
    */

    tag "$sampleID"

    conda params.mtbseq_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_mtbseq
            else if (workflow.containerEngine == 'docker') return params.docker_mtbseq
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_mtbseq
            else return null
        }
    
    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}", mode: 'copy'

    input:
        tuple val(sampleID), 
                path(mtbc_forward), //, stageAs: "${sampleID}_R1.fastq.gz"
                path(mtbc_reverse), //, stageAs: "${sampleID}_R2.fastq.gz"
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
        path("Mpileup/*")
        path("Position_Tables/*")
        path("Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab")

        // tuple for updating the sample_ch
        tuple val(sampleID), 
                path(mtbc_ont_reads), 
                path("Classification/${sampleID}.Strain_Classification.tab"), 
                path("Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab"), 
                path("Position_Tables/${sampleID}.gatk_position_table.tab"), 
                path("Called/${sampleID}.gatk_position_variants_*.tab"),  
                path(tbdb_out), path(who_out), 
                path(mtbseq_vcf),path("Mpileup/${sampleID}.gatk.mpileup"),              emit: updated_sample_ch4

    script:

        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        # Rename the reads to the intended name without the "_mtbc"
            mtbc_Fname1=\$(ls ${mtbc_ont_reads})
            mtbc_Fname2=\$(ls ${mtbc_ont_reads} | sed 's/_mtbc//g')
            mv \${mtbc_Fname1} \${mtbc_Fname2}

        # map reads to reference
            minimap2 -ax-ont ${mtbc_ont_reads} 

        # Run MTBseq for a single sample

            MTBseq --step TBrefine --continue \\
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

        # restore the symbolic link names
            mv \${mtbc_Fname2} \${mtbc_Fname1}; mv \${mtbc_Rname2} \${mtbc_Rname1}

        # Rename the stats and class outputs to have unique names
        ## this prevent clashes later on
            cp Classification/Strain_Classification.tab Classification/${sampleID}.Strain_Classification.tab
            cp Statistics/Mapping_and_Variant_Statistics.tab Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab

        # remove the published reads from the previous module:
            rm -f  ${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}_mtbc_R1.fastq.gz       
            rm -f  ${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}_mtbc_R2.fastq.gz
            rm -f  ${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}/tbdb-${sampleID}.results.txt
        """

}