process MTBSEQ_SINGLE {

    tag "$sampleID"

    conda params.mtbseq_env

    container { if (workflow.containerEngine == 'singularity') { 
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
            } else { 'quay.io/biocontainers/mtbseq' }
    }
    
    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}", mode: 'copy'

    input:
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
                path(tbdb_out), path(who_out), path(mtbseq_vcf)
                
    output:
        path("Bam/*")
        path("Called/*")
        path("Classification/${sampleID}.Strain_Classification.tab")
        path("GATK_Bam/*")
        path("Mpileup/*")
        path("Position_Tables/*")
        path("Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab")

        // tuple for updating the sample_ch
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse), 
                path("Classification/${sampleID}.Strain_Classification.tab"), 
                path("Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab"), 
                path("Position_Tables/${sampleID}.gatk_position_table.tab"), 
                path("Called/${sampleID}.gatk_position_variants_*.tab"),  
                path(tbdb_out), path(who_out), 
                path(mtbseq_vcf),path("Mpileup/${sampleID}.gatk.mpileup"),              emit: updated_sample_ch4

    script:

        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        # Brfiedly rename the reads to the intended name without the "_mtbc"
            mtbc_Fname1=\$(ls ${mtbc_forward}); mtbc_Rname1=\$(ls ${mtbc_reverse})
            mtbc_Fname2=\$(ls ${mtbc_forward} | sed 's/_mtbc//g'); mtbc_Rname2=\$(ls ${mtbc_reverse} | sed 's/_mtbc//g')
            mv \${mtbc_Fname1} \${mtbc_Fname2}; mv \${mtbc_Rname1} \${mtbc_Rname2}


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

        # restore the symbolic link names
            mv \${mtbc_Fname2} \${mtbc_Fname1}; mv \${mtbc_Rname2} \${mtbc_Rname1}

        # Rename the stats and class outputs to have unique names
        ## this prevent clashes later on
            mv Classification/Strain_Classification.tab Classification/${sampleID}.Strain_Classification.tab
            mv Statistics/Mapping_and_Variant_Statistics.tab Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab

        """

}