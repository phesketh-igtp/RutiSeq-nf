process MTBSEQ_LINEAGE_PAIRWISE {

    tag " ${runID}: ${lineage} "

    conda params.mtbseq_env

    // TODO: container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
    ///        } else { 'quay.io/biocontainers/mtbseq' }
    ///}
                
    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/", mode: 'copy'

    input:
        val runID
        tuple val(lineage), val(sampleIDs)

    output:
        path("${lineage}_samples.txt")

        // Amend outputs
        tuple val(lineage), path("Amend/"),                                             emit: amend_dir
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.fasta"),
                            path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.tab"),   emit: snp_phylogeny_ch
        
        path("Amend/")                    
        path("Amend/${lineage}_joint_*_amended.tab")
        path("Amend/${lineage}_joint_*_amended_u*_phylo.fasta")
        path("Amend/${lineage}_joint_*_amended_u*_phylo.plainIDs.fasta")
        path("Amend/${lineage}_joint_*_amended_u*_phylo.tab")
        path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.fasta")
        path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.plainIDs.fasta")
        path("Amend/${lineage}_joint_*_amended_u*_phylo_w*_removed.tab")
        path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.tab")

        // Join output
        tuple val(lineage), path("Join/"),                                               emit: join_dir                                             
        path("Join/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.log") 
        path("Join/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.tab")

        // Groups
        tuple val(lineage), path("Group/"),                                              emit: group_dir
        path("Groups/${lineage}_joint_*.groups")
        path("Groups/${lineage}_joint_*.matrix")
        path("Groups/${lineage}_*.list")
        path("Groups/${lineage}_clusters.tsv"),                                         emit: clusters

        //Matrix ouput

        path("Matrices/${lineage}.matrix"),                                              emit: matrix_dir

    script:

    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    # make the expected directories
        mkdir Position_Tables/ Called/ Groups/ Matrices

    # create the list of the sampleIDs within that lineage
        echo "${sampleIDs.join('\n')}" > samplesID.list

        sed 's@_@\t@g' samplesID.list > ${lineage}_samples.txt

    # Create symbolic links to the apprpriate files
        while IFS=',' read -r samples; do

            for file in ${params.outdir}/bbdd/mtbseq/samples/\${samples}/Position_Tables/*.tab; do
                ln -s \${file} Position_Tables/\$(basename \$file); done

            for file in ${params.outdir}/bbdd/mtbseq/samples/\${samples}/Called/*.tab; do
                ln -s \${file} Called/\$(basename \$file); done

        done < samplesID.list

        MTBseq --step TBjoin \\
            --thread ${task.cpus} \\
            --project ${lineage} \\
            --samples samplesID.list \\
            ${additional_args} \\
            1>>.command.out \\
            2>>.command.err \\
            || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

        MTBseq --step TBamend \\
            --thread ${task.cpus} \\
            --project ${lineage} \\
            --samples samplesID.list \\
            ${additional_args} \\
            1>>.command.out \\
            2>>.command.err \\
            || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    # Get the list of SNP distances to analyse
        echo '${params.mtbseq_snp_distance.join('\n')}' > snp_distances

        while read -r distance; do

                MTBseq --step TBgroups \\
                    --thread ${task.cpus} \\
                    --project ${lineage} \\
                    --samples samplesID.list \\
                    --distance \${distance} \\
                    ${additional_args} \\
                    1>>.command.out \\
                    2>>.command.err \\
                    || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

            # Wrangle the group list into a useful format
                sed '0,/### Output as lists:/d' Groups/${lineage}_joint_*_\${distance}.groups \\
                                                            > Groups/${lineage}_\${distance}.list
                
                # Add lineage and SNP distance the tsv file
                sed -i "s@^@${lineage}\t\${distance}\t@g" Groups/${lineage}_\${distance}.list

        done < snp_distances

    # Move and rename the matrices
        mv Groups/${lineage}*.matrix Matrices/${lineage}.matrix

        cat Groups/${lineage}_*.list > Groups/${lineage}_clusters.tsv


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

        cat Matrices/${lineage}.matrix.head Matrices/${lineage}.matrix > Matrices/${lineage}.matrix.tsv

        # remove the intermediates
        rm Matrices/${lineage}.matrix.head Matrices/${lineage}.matrix

    """
}