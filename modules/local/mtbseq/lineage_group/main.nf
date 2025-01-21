process MTBSEQ_LINEAGE_GROUP {

    tag "${runID}: ${lineage} | d${distance}"

    conda params.mtbseq_env

    // TODO: container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
    ///        } else { 'quay.io/biocontainers/mtbseq' }
    ///}
                
    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/", mode: 'copy'

    input:
        val runID
        tuple val(lineage), val(distance), 
                path(joint_dir), 
                path(amend_dir),
                path(sample_txt)

    output:
        // Groups
        path("Group/*")
        path("Groups/${lineage}_d${distance}.groups")
        path("Groups/${lineage}_clusters.d${distance}.tsv"),            emit: clusters

        //Matrix ouput
        path("Matrices/*")
        path("Matrices/${lineage}.d${distance}.matrix.tsv"),            emit: matrix_dir

    script:

        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        mkdir -p Groups/

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
                ${additional_args} \\
                1>>.command.out \\
                2>>.command.err \\
                || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

        # Rename the groups file for simplicity
            cat Groups/${lineage}_*d${distance}.groups > Groups/${lineage}_d${distance}.groups

        # Wrangle the group list into a useful format
            cat Groups/${lineage}_joint_*_${distance}.groups \\
                    | sed '0,/### Output as lists:/d' \\
                        > Groups/${lineage}_d${distance}.clusters.tsv
                    
        # Add lineage and SNP distance the tsv file
            sed -i "s@^@${lineage}\t${distance}\t@g" Groups/${lineage}_${distance}.list

        # Move and rename the matrices
            mv Groups/${lineage}*.matrix Matrices/${lineage}.d${distance}.matrix

            cat Groups/${lineage}_*.list > Groups/${lineage}_clusters.d${distance}.tsv

        ## Correct the format of the matrix for importing to R

            # cut the headers column from the matrix
            cut -f1 Matrices/${lineage}.matrix > Matrices/${lineage}.matrix

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
                ' Matrices/${lineage}.matrix.ids | sed 's/^/sampleID\t/g' > Matrices/${lineage}.matrix.head

                cat Matrices/${lineage}.matrix.head Matrices/${lineage}.matrix > Matrices/${lineage}.d${distance}.matrix.tsv

            # remove the intermediates
                rm Matrices/${lineage}.matrix.head Matrices/${lineage}.matrix

                gzip --best Matrices/${lineage}.d${distance}.matrix.tsv

        """
}