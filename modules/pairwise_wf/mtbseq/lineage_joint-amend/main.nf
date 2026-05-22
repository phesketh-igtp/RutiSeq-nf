process MTBSEQ_LINEAGE_JOINT_AMEND {

    tag "${lineage} "

    conda params.mtbseq_env

    // TODO: container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
    ///        } else { 'quay.io/biocontainers/mtbseq' }
    ///}
                
    publishDir "${params.outDir}/db/comparison/mtbseq/${lineage}/",         
        mode: 'copy', 
        overwrite: true

    input:
        tuple val(lineage), 
            val(sampleIDs),
            val(sampleID_count)

    output:
    // Main Outputs
        path("${lineage}_samples.txt")
        path("Amend/*") 
        path("Joint/*")   

    // snp_phylogeny_ch
        tuple val(lineage),
            path("Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta"),
            path("Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo.tab"),   
            emit: snp_phylogeny_ch

    // Tuple for creating MTBSeq --step TBgroup channel
        path("mtbseq-group.tuple.csv"), emit: mtbseq_group_tuple_csv
    
    script:

    // Defined variables
        def additional_args = task.ext.additional_args ?: ''
        def OUTDIR = "${params.outDir}/db/comparison/mtbseq/${lineage}"
        def MTBSEQ_FASTA_OUT = "Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta"
        def MTBSEQ_TAB_OUT = "Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}_amended_u${params.mtbseq_unambig}_phylo.tab"
        def snp_distances = params.mtbseq_snp_distance.join('\n')

        """
        # make the expected directories
            mkdir -p Position_Tables/ Called/ Amend/ Joint/

        # create the list of the sampleIDs within that lineage
            echo "${sampleIDs.join('\n')}" | sort | uniq > samplesID.list
            sed 's@_@\t@g' samplesID.list > ${lineage}_samples.txt

        # Get list of sample names from the old version:
            grep '>' ${OUTDIR}/${MTBSEQ_FASTA_OUT} | sed 's@>@@g' > ${lineage}_samples.old.txt

        # Check the sampleIDs all match to decide to run or skip TBjoin/TBamend
        if cmp -s "samplesID.list" "${lineage}_samples.old.txt"; then

            echo "Skipping MTBseq Join and Amend steps for lineage: ${lineage}. Linking .
                    ${params.outDir}/db/comparison/mtbseq/${lineage}/Joint/ and 
                    ${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/"

            # Should be:
            ln -s ${params.outDir}/db/comparison/mtbseq/${lineage}/Joint/* Joint/
            ln -s ${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/* Amend/

        else

            echo "Running MTBseq Join and Amend steps for lineage: ${lineage}"

            # Remove any previously generated files to generate new ones
                rm -rf ${OUTDIR}/Joint/* \\
                        ${OUTDIR}/Amend/* \\
                        ${OUTDIR}/Groups/*

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

            ## Check that there are AT MINIMUM 4 MILLION positions in the Joint table
            if [[ $(cut -f1 "$MTBSEQ_TAB_OUT" | tail -1) -lt 4000000 ]]; then
                echo -e "$MTBSEQ_TAB_OUT is likely truncated. Revisit individual samples"
                exit 1
            fi

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
        fi

        ### create lineage csv for creating new channels
        echo '${snp_distances}' > snp_distances
        for distance in \$(cat snp_distances); do
            echo "${lineage},\${distance},${params.outDir}/db/comparison/mtbseq/${lineage}/Joint,${params.outDir}/db/comparison/mtbseq/${lineage}/Amend,${params.outDir}/db/comparison/mtbseq/${lineage}/${lineage}_samples.txt,${sampleID_count}" >> mtbseq-group.tuple.csv
        done

        # Check the Amend file is it actually has any sequences in the fasta
        sum_len=\$(seqkit stats -T ${MTBSEQ_FASTA_OUT} | sed '1d' | head -1 | cut -f5)
        if [[ \${sum_len} == 0 ]]; then
            echo "Error: The Amend/ FASTA file ${MTBSEQ_FASTA_OUT} is empty." >&2
            exit 1
        fi
        """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: v2.1.0
@description:
    This process runs the MTBseq Join and Amend steps for a given lineage.
    It takes the lineage and sample IDs as input and produces the output files
    for the MTBseq Join and Amend steps.
    It also creates a lineage_samples.txt file for the MTBseq Join step.
@changelog:
    v1.0.0-2025-04-01: Initial version
    v1.1.0-2025-04-09: Changed - Removed the zipping of the files (caused issues)
    v2.0.0-2026-04-24: Added full parameters to outputs names to ensure accuracy. Included the count of sampleIDs in the channel.
    v2.0.1-2026-04-27: Moved snp distances to the definition line before the shell script
    v2.1.0-2026-05-22: Changed to IF statement to be a comparison of the DEFLINES from the old and new files
                        This bypasses cases when the number of sample remains the same but different samples are in the 
                        collection (remove + replace 1)
*/
