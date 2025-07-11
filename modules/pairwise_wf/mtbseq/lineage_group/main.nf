process MTBSEQ_LINEAGE_GROUP {

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
*/

    tag "${lineage}; t=${distance}"

    conda params.mtbseq_env

    // TODO: container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
    ///        } else { 'quay.io/biocontainers/mtbseq' }
    ///}
                
    publishDir "${params.outDir}/bbdd/mtbseq/pairwise/${lineage}/", mode: 'copy'

    input:
        val runID
        tuple val(lineage), 
                val(distance), 
                path(joint_dir), 
                path(amend_dir),
                path(sample_txt)

    output:
        // Groups
        path("Groups/*")
        path("Groups/${lineage}_d${distance}.clusters.tsv"), emit: clusters

        //Matrix ouput
        path("Matrices/*")
        path("Matrices/${lineage}.d${distance}.matrix.tsv"), emit: matrix_dir

        path("${lineage}_handover.txt"),                     emit: handover

    script:

        //def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        mkdir -p Groups/ Matrices

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

        # Rename the groups file for simplicity
            cat Groups/${lineage}_*_d${distance}.groups > Groups/${lineage}_d${distance}.groups

        # Wrangle the group list into a useful format
            cat Groups/${lineage}_d${distance}.groups \\
                    | sed '0,/### Output as lists:/d' \\
                        > Groups/${lineage}_d${distance}.clusters.tsv
                    
        # Add lineage and SNP distance the tsv file
            sed -i "s@^@${lineage}\t${distance}\t@g" Groups/${lineage}_d${distance}.clusters.tsv

        # Move and rename the matrices
            mv Groups/${lineage}*.matrix Matrices/${lineage}.d${distance}.matrix

        ## Correct the format of the matrix for importing to R

            # cut the headers column from the matrix
            cut -f1 Matrices/${lineage}.d${distance}.matrix > Matrices/tmp.${lineage}.matrix.ids

            # transpose the first colum long to wide (tab seperated)
                awk '                                         
                {
                    for (i = 1; i <= NF; i++) {
                        arr[i] = (arr[i] ? arr[i] "\t" : "") \$i;
                    }
                }
                END {
                    for (i = 1; i in arr; i++) {
                        print arr[i];
                    }
                }
                ' Matrices/tmp.${lineage}.matrix.ids | sed 's/^/sampleID\t/g' > Matrices/tmp.${lineage}.matrix.head

                cat Matrices/tmp.${lineage}.matrix.head Matrices/${lineage}.d${distance}.matrix > Matrices/${lineage}.d${distance}.matrix.tsv

            # remove the intermediates
                rm Matrices/tmp.${lineage}.*

        # Handover to prevent reruns
            echo "${lineage} handover" > ${lineage}_handover.txt

        """
}