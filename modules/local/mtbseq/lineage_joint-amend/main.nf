process MTBSEQ_LINEAGE_JOINT_AMEND {

    tag " ${runID}: ${lineage} "

    conda params.mtbseq_env

    // TODO: container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
    ///        } else { 'quay.io/biocontainers/mtbseq' }
    ///}
                
    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/", mode: 'copy', overwrite: true

    input:
        val runID
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
        path("Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo.tab")

        // Join output
        tuple val(lineage), path("Joint/*")      
        path("Joint/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.log") 
        path("Joint/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.tab")
        
        // Tuple for creating MTBSeq --step TBgroup channel
        path("mtbseq-group.tuple.csv"),     emit: mtbseq_group_tuple_csv
    
    script:

        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        rm -rf ${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Amend/*
        gzip --quiet --best ${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Joint/*

        # make the expected directories
            mkdir -p Position_Tables/ Called/ Amend/ Joint/

        # create the list of the sampleIDs within that lineage
            echo "${sampleIDs.join('\n')}" | sort | uniq > samplesID.list
            sed 's@_@\t@g' samplesID.list > ${lineage}_samples.txt

        # Create symbolic links to the apprpriate files (only if the file does not exist)
            while IFS=',' read -r samples; do

                for file in ${params.outdir}/bbdd/mtbseq/samples/\${samples}/Position_Tables/*.tab; do
                    dest="Position_Tables/\$(basename "\$file")"
                    if [[ ! -e "\$dest" && ! -L "\$dest" ]]; then
                        ln -s "\$file" "\$dest"
                    fi
                done

                for file in ${params.outdir}/bbdd/mtbseq/samples/\${samples}/Called/*.tab; do
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

        ### create lineage csv for creating new channels
            echo '${params.mtbseq_snp_distance.join('\n')}' > snp_distances

            for distance in \$(cat snp_distances); do
                echo "${lineage},\${distance},${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Joint,${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Amend,${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/${lineage}_samples.txt" >> mtbseq-group.tuple.csv
            done                
        """
}