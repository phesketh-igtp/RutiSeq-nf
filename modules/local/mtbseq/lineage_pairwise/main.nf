process MTBSEQ_LINEAGE_PAIRWISE {

    tag "${runID}: ${lineage}"

    conda params.mtbseq_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
            } else { 'quay.io/biocontainers/mtbseq' }
    }
                
    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/${lineage}", mode: 'copy'

    input:
        val runID
        tuple val(lineage), val(sampleIDs)

    output:
        // Amend outputs
        path("${lineage}_samples.txt")
        tuple val(lineage), path("Amend/"),                                     emit: amend_dir                                   
        path("Amend/${lineage}_joint_*_amended.tab")
        path("Amend/${lineage}_joint_*_amended_u*_phylo.fasta")
        path("Amend/${lineage}_joint_*_amended_u*_phylo.plainIDs.fasta")
        path("Amend/${lineage}_joint_*_amended_u*_phylo.tab")
        path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.fasta")
        path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.plainIDs.fasta")
        path("Amend/${lineage}_joint_*_amended_u*_phylo_w*_removed.tab")
        path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.tab")
        // Join output
        tuple val(lineage), path("Join/"),                                      emit: join_dir                                             
        path("Join/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.log") 
        path("Join/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.tab")                   

    script:

    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    # make the expected directories
        mkdir Position_Tables/ Called/

    # create the list of the sampleIDs within that lineage
        echo "${sampleIDs.join('\n')}" > samplesID.list

        sed 's@\t@_@g' samplesID.list > ${lineage}_samples.txt

        while IFS=',' read -r sampleID; do
            ln -s ${params.outdir}/bbdd/mtbseq/pairwise/\${sampleID}/Position_Tables/* Position_Tables/
            ln -s ${params.outdir}/bbdd/mtbseq/pairwise/\${sampleID}/Called/* Called/
        done < samplesID.list

        MTBseq --step TBjoin \\
            --thread ${task.cpus} \\
            --project ${lineage} \\
            --sample samplesID.list \\
            ${additional_args} \\
            1>>.command.out \\
            2>>.command.err \\
            || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

        MTBseq --step TBamend \\
            --thread ${task.cpus} \\
            --project ${lineage} \\
            --sample samplesID.list \\
            ${additional_args} \\
            1>>.command.out \\
            2>>.command.err \\
            || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    # Get the list of SNP distances to analyse
        echo '${params.lineage_pairwise.join('\n')}' > snp_distances

    while read -r distance; do

            MTBseq --step TBgroup \\
                --thread ${task.cpus} \\
                --project ${lineage} \\
                --sample samplesID.list \\
                --distance \${distance} \\
                ${additional_args} \\
                1>>.command.out \\
                2>>.command.err \\
                || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    done < snp_distances

    # Modify the matrix file to be usable for future:
    

    """
}