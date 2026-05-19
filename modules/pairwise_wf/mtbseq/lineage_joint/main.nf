process MTBSEQ_LINEAGE_JOINT {

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
    
        tuple val(lineage), 
            val(sampleIDs),
            val(sampleID_count),
            path("Joint/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}.tab"),
            path("Joint/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples${sampleID_count}.log"), 
            emit: mtbseq_joint_out
    
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

        if [[ ! -f "${OUTDIR}/${MTBSEQ_FASTA_OUT}" && ! -f "${OUTDIR}/${MTBSEQ_TAB_OUT}" ]]; then

            ## The correct Amend/ file does not exist, so we need to run MTBseq Join and Amend steps
                echo "Running MTBseq Join and Amend steps for lineage: ${lineage}"
                
            # Remove any previously generated files to generate new ones
                rm -rf ${OUTDIR}/Joint/*

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

            # Join the MTBseq Position Tables using custom Python script
            python ${params.scriptDir}}/py/pl.TBjoin.py \
                --var-dir /imppc/labs/emlab/share/GitHub/RutiSeq-nf/db/H37Rv_db/ \
                --pos-out Position_Tables \
                --call-out Called \
                --join-out Joint \
                --group-name ${lineage} \\
                --ref M._tuberculosis_H37Rv_2015-11-13.fasta \\
                --refg M._tuberculosis_H37Rv_2015-11-13_genes.txt \\
                --mincovf ${params.mtbseq_mincovf} \\
                --mincovr ${params.mtbseq_mincovr} \\
                --minphred20 ${params.mtbseq_minphred20}  \\
                --minfreq ${params.mtbseq_minfreq}

        else

            echo "Skipping MTBseq --step TBjoin steps for lineage: ${lineage}. Delete the existing files to re-run.
                    ${params.outDir}/db/comparison/mtbseq/${lineage}/Joint/"

            # Should be:
            cp ${params.outDir}/db/comparison/mtbseq/${lineage}/Joint/* Joint/
        
        fi
        """
}

/*
@author: Poppy J Hesketh Best
@date: 2026-05-19
@version: v1.0.0
@description:
    This process runs the MTBseq Join and Amend steps for a given lineage.
    It takes the lineage and sample IDs as input and produces the output files
    for the MTBseq Join and Amend steps.
    It also creates a lineage_samples.txt file for the MTBseq Join step.
@changelog:
    v1.0.0-2026-05-19: Initial version

TODO: Would be better to compare the sampleIDs from the fasta files and 
            if they are differente, then remove and repeat the analysis, or if they are the same
            continue
*/
