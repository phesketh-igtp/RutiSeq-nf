process MTBSEQ_LINEAGE_JOINT_AMEND {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process runs the MTBseq Join and Amend steps for a given lineage.
        It takes the lineage and sample IDs as input and produces the output files
        for the MTBseq Join and Amend steps.
        It also creates a lineage_samples.txt file for the MTBseq Join step.
    @changelog:
        v1.0.0-2025-04-01: Initial version
        v1.1.0-2025-04-09: Changed - Removed the zipping of the files (caused issues)

    TODO: Would be better to compare the sampleIDs from the fasta files and 
            if they are differente, then remove and repeat the analysis, or if they are the same
            continue
*/

    tag "${params.runID}: ${lineage} "

    conda params.mtbseq_env

    // TODO: container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
    ///        } else { 'quay.io/biocontainers/mtbseq' }
    ///}
                
    publishDir "${params.outDir}/db/comparison/mtbseq/${lineage}/",         
        mode: 'copy', 
        overwrite: true

    input:
        tuple val(lineage), val(sampleIDs)

    output:
        path("${lineage}_samples.txt")

        // Amend outputs
        tuple val(lineage), path("Amend/*")
        tuple val(lineage), path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta"),
                            path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo.tab"),   
                                            emit: snp_phylogeny_ch

        path("Amend/*")                    
        path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended.tab")
        path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo.fasta")
        path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo.plainIDs.fasta")
        path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta")
        path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.plainIDs.fasta")
        path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}_removed.tab")
        path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.tab")

        // Join output
        tuple val(lineage), path("Joint/*")      
        path("Joint/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.log") 
        path("Joint/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.tab")
        
        // Tuple for creating MTBSeq --step TBgroup channel
        path("mtbseq-group.tuple.csv"),     emit: mtbseq_group_tuple_csv
    
    script:

    def additional_args = task.ext.additional_args ?: ''
    def sampleID_count = sampleIDs.size() // Get the count of sample IDs
    def mtbseq_fasta_out = "Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta"
    
        """
        # make the expected directories
            mkdir -p Position_Tables/ Called/ Amend/ Joint/

        # create the list of the sampleIDs within that lineage
            echo "${sampleIDs.join('\n')}" | sort | uniq > samplesID.list
            sed 's@_@\t@g' samplesID.list > ${lineage}_samples.txt

        if [[ ! -f "${params.outDir}/db/comparison/mtbseq/${lineage}/${mtbseq_fasta_out}" ]]; then

            ## The correct Amend/ file does not exist, so we need to run MTBseq Join and Amend steps
                echo "Running MTBseq Join and Amend steps for lineage: ${lineage}"
                
            # Remove any previously generated files to generate new ones
                rm -rf ${params.outDir}/db/comparison/mtbseq/${lineage}/Joint/*
                rm -rf ${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/*
                rm -rf ${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/*
                rm -rf ${params.outDir}/db/comparison/mtbseq/${lineage}/Groups/*
                rm -rf ${params.outDir}/db/comparison/mtbseq/${lineage}/Matrices/*

            # Create symbolic links to the apprpriate files (only if the file does not exist)
                    while IFS=',' read -r samples; do

                        for file in ${params.outDir}/db/samples/\${samples}/mtbseq/Position_Tables/*.tab; do
                            dest="Position_Tables/\$(basename "\$file")"
                            if [[ ! -e "\$dest" && ! -L "\$dest" ]]; then
                                ln -s "\$file" "\$dest"
                            fi
                        done

                        for file in ${params.outDir}/db/samples/\${samples}/mtbseq/Called/*.tab; do
                            dest="Called/\$(basename "\$file")"
                            if [[ ! -e "\$dest" && ! -L "\$dest" ]]; then
                                ln -s "\$file" "\$dest"
                            fi
                        done

                    done < samplesID.list

            ## MTBseq Join using the first SNP distance
                MTBseq --step TBjoin ${additional_args} \\
                        --thread        ${task.cpus} \\
                        --project       ${lineage} \\
                        --samples       ${lineage}_samples.txt \\
                        --minbqual      ${params.mtbseq_minbqual} \\
                        --mincovf       ${params.mtbseq_mincovf} \\
                        --mincovr       ${params.mtbseq_mincovr} \\
                        --minphred20    ${params.mtbseq_minphred20} \\
                        --minfreq       ${params.mtbseq_minfreq} \\
                        --unambig       ${params.mtbseq_unambig} \\
                        --window        ${params.mtbseq_window} \\
                            1>>.command.out \\
                            2>>.command.err || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

            ## MTBseq Join using the first SNP distance
                    MTBseq --step TBamend ${additional_args} \\
                        --thread        ${task.cpus} \\
                        --project       ${lineage} \\
                        --samples       ${lineage}_samples.txt \\
                        --minbqual      ${params.mtbseq_minbqual} \\
                        --mincovf       ${params.mtbseq_mincovf} \\
                        --mincovr       ${params.mtbseq_mincovr} \\
                        --minphred20    ${params.mtbseq_minphred20} \\
                        --minfreq       ${params.mtbseq_minfreq} \\
                        --unambig       ${params.mtbseq_unambig} \\
                        --window        ${params.mtbseq_window} \\
                            1>>.command.out \\
                            2>>.command.err || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq       

        else

            echo "Skipping MTBseq Join and Amend steps for lineage: ${lineage}. Delete the existing files to re-run.
                    ${params.outDir}/db/comparison/mtbseq/${lineage}/Joint/ and 
                    ${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/"

            # Should be:
            cp ${params.outDir}/db/comparison/mtbseq/${lineage}/Joint/* Joint/
            cp ${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/* Amend/
        
        fi

        ### create lineage csv for creating new channels
                echo '${params.mtbseq_snp_distance.join('\n')}' > snp_distances

                for distance in \$(cat snp_distances); do
                    echo "${lineage},\${distance},${params.outDir}/db/comparison/mtbseq/${lineage}/Joint,${params.outDir}/db/comparison/mtbseq/${lineage}/Amend,${params.outDir}/db/comparison/mtbseq/${lineage}/${lineage}_samples.txt" >> mtbseq-group.tuple.csv
                done

        # Check the Amend file is it actually has any sequences in the fasta
            sum_len=\$(seqkit stats -T ${mtbseq_fasta_out} | sed '1d' | head -1 | cut -f5)
                if [[ \${sum_len} == 0 ]]; then
                    echo "Error: The Amend/ FASTA file ${mtbseq_fasta_out} is empty." >&2
                    exit 1
                fi
        """
}

/*
# Define the fasta output
    #sampleID_count=\$(wc -l samplesID.list | cut -f1 -d' ')
    #mtbseq_fasta_out=\$(echo -e "Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta")

*/