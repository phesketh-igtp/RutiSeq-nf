process MTBSEQ_LINEAGE_GROUP {

    tag "${lineage}; t=${distance}"

    conda params.mtbseq_env

    // TODO: container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
    ///        } else { 'quay.io/biocontainers/mtbseq' }
    ///}
                
    publishDir "${params.outDir}/db/comparison/mtbseq/${lineage}/", 
        mode: 'copy', 
        overwrite: true
    //storeDir "${params.outDir}/db/comparison/mtbseq/${lineage}/"

    input:
        tuple val(lineage), 
                val(distance), 
                path(joint_dir), 
                path(amend_dir),
                path(sample_txt),
                val(sampleID_count)

    output:
        //Matrix ouput
        path("Groups/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}_d${distance}.groups")
        path("Groups/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.matrix")

        // Nexus output
        tuple val(lineage), 
            val(distance), 
            path("${params.runID}_${lineage}_snps.fasta"),
            path("${params.runID}_${lineage}_snps.tab"),
            path("Groups/${params.runID}_${lineage}_d${distance}.clusters.tsv"), 
            path("Groups/${params.runID}_${lineage}_d${distance}.singletons.tsv"), emit: clusters

    script:

        //def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file
        //Groups/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window].matrix
        def groups_tab = "Groups/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}_d${distance}.groups"
        def groups_mat = "Groups/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.matrix"
        def amend_fasta = "Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta"
        def amend_tab = "Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.tab"

        """
        # Make output directories (local)
        mkdir -p Groups/ Matrix/

        if [[ -f ${params.outDir}/db/comparison/mtbseq/${lineage}/${groups_tab} && \
            -f ${params.outDir}/db/comparison/mtbseq/${lineage}/${groups_mat} ]]; then

            ln -s ${params.outDir}/db/comparison/mtbseq/${lineage}/${groups_tab} Groups/
            ln -s ${params.outDir}/db/comparison/mtbseq/${lineage}/${groups_mat} Groups/

        else

            # Clean up the last results
            rm -rf ${params.outDir}/db/comparison/mtbseq/${lineage}/Groups/*

            ## MTBseq TBgroups using the first SNP distance
            MTBseq --step TBgroups \\
                --thread        ${task.cpus} \\
                --project       ${lineage} \\
                --samples       ${sample_txt} \\
                --distance      ${distance} \\
                --minbqual      ${params.mtbseq_minbqual} \\
                --mincovf       ${params.mtbseq_mincovf} \\
                --mincovr       ${params.mtbseq_mincovr} \\
                --minphred20    ${params.mtbseq_minphred20} \\
                --minfreq       ${params.mtbseq_minfreq} \\
                --unambig       ${params.mtbseq_unambig} \\
                --window        ${params.mtbseq_window} \\
                    1>>.command.out \\
                    2>>.command.err \\
                    || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq
        fi

        # Wrangle the group list into a useful format
        cat ${groups_tab} \\
                | sed '0,/### Output as lists:/d' \\
                | sed "s@^@${lineage}\tdist_${distance}\t@g" \\
                    > tmp.${lineage}_d${distance}.clusters.tsv

        grep 'group_' tmp.${lineage}_d${distance}.clusters.tsv \\
                | sed "s@group_@@g" \\
                | sed s/\$/-${distance}-${lineage}/g \\
                | sed 's/-lineage/-L/g' \\
                > Groups/${params.runID}_${lineage}_d${distance}.clusters.tsv

        grep 'ungrouped' tmp.${lineage}_d${distance}.clusters.tsv \\
                | sed "s@ungrouped@singleton@g" > Groups/${params.runID}_${lineage}_d${distance}.singletons.tsv

        cat ${amend_fasta} > ${params.runID}_${lineage}_snps.fasta
        cat ${amend_tab} > ${params.runID}_${lineage}_snps.tab
        """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 1.0.1
@description:
    This process runs the MTBseq TBgroups step on the joint and amend directories
        for each lineage and distance. It takes the output from the
        MTBSEQ_LINEAGE_JOINT_AMEND() process and runs the TBgroups step on the joint
        and amend directories. It also renames the output files for simplicity.
        It also wrangles the output matrix into a useful format for downstream analysis.
        Occasionally, the MTBseq TBgroups step will fail and produice empty files/no-files.
@last_updated: 2025-04-01
@changelog:
    v1.0.0-2025-04-01: Initial version + documnetadocumentation
    v1.0.1-2025-05-19: Removed additonal_args as it was not utilised and was creating inconsistencies
    v1.0.2-2026-05-19: Removed the matrix output.
*/