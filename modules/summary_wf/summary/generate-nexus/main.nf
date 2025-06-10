process GENERATE_NEXUS {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process generates a NEXUS file for each cluster of genomes
        using the SNPs from the reference genome and the MTB_anc
        (Valencian ancestor, Iñaki Comas 2013) as a reference. The resulting 
        NEXUS file is intended to be use with PopArt for visualising median-joining
        networks.
*/

    conda params.snp_profiling_env 

    tag "cluster: ${clusterID}"

    publishDir "${params.outDir}/results/${runID}/networks/", mode: 'copy'

    input:
        val(runID)
        file(pairwise_clusters)
        tuple val(lineage), 
                val(clusterID),
                file(snp_fasta),
                file(snp_tab)

    output:
        path("fasta/*"),        optional: true
        path("positions/*"),    optional: true
        path("nexus/*"),        optional: true

        tuple val(lineage), 
                val(clusterID), 
                path("fasta/${clusterID}_refseq.fasta"),
                path("positions/${clusterID}_genomic_positions.tab"),
                path("${clusterID}.snp.tab"),   emit: variant_sites_for_tabulation

    script:

        """
        # run the nexus script
            bash ${params.script_dir}/shell/create-variable-region-nexus.sh \\
                    -c ${clusterID} \\
                    -p ${pairwise_clusters} \\
                    -f ${snp_fasta} \\
                    -t ${snp_tab} \\
                    -m ${params.mtbc_ancestor_path} \\
                    1>>.command.out \\
                    2>>.command.err || true

        # simplify the name of the file
            cp ${snp_tab} ${clusterID}.snp.tab

        # check if the nexus generation was successful
            if [[ ! -f fasta/${clusterID}_refseq.fasta ]]; then echo "Nexus generation failed" > fasta/${clusterID}_refseq.fasta; fi
            if [[ ! -f positions/${clusterID}_genomic_positions.tab ]]; then echo "Nexus generation failed" > positions/${clusterID}_genomic_positions.tab; fi

        """

}