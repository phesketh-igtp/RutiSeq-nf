process MTBSEQ_LINEAGE_GROUPING {

    tag "${runID}_${lineage}"

    conda params.mtbseq_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
            } else { 'quay.io/biocontainers/mtbseq' }
    }

    publishDir "${params.outdir}/bbdd/mtbseq/lineages/${lineage}", mode: 'copy'

    input:
        val runID
        tuple val(lineage), path(samples_file)
        path amend_dir
        path join_dir


    output:
        path("Groups/")
        path("Groups/${lineage}_joint_*")
        path("Groups/${lineage}_clusters.tsv"),     emit: clusters
        path("Matrices/${lineage}_d*.matrix"),      emit: matrices
        

    script:

        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        mkdir -p Groups/ Matrices/

        # Create sample list of all the MTB lineages to be analyzed

        for snp_distance in ${params.snp_distance.join(' ')}; do

            MTBseq --step TBgroup \\
                --thread ${task.cpus} \\
                --prefix ${lineage} \\
                --sample ${samples_file} \\
                ${additional_args} \\
                --distance \${snp_distance} \\
                1>>.command.out \\
                2>>.command.err \\
                || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

            # Create the tidy list of genomes belonging to a cluster
            cat Groups/${lineage}_joint_*d\${snp_distance}.groups > Groups/${lineage}_d\${snp_distance}.list
            sed -i '0,/### Output as lists:/d' Groups/${lineage}_d\${snp_distance}.list
            sed -i "s@^@${lineage}\t\${snp_distance}\t@g" Groups/${lineage}_d\${snp_distance}.list

            # Rename the matrix file to have the distance value
            mv Groups/${lineage}_joint_*.matrix Matrices/${lineage}_d\${snp_distance}.matrix

        done

        # Combine all distance-specific lists into a single clusters file
        cat Groups/${lineage}_d*.list > Groups/${lineage}_clusters.tsv

        """
}