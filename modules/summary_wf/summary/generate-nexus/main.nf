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

    conda params.snippy_env 

    tag "cluster: ${clusterID}"

    publishDir "${params.outDir}/results/${params.runID}/networks/", mode: 'copy'

    input:
        tuple val(lineage), 
                val(clusterID),
                file(snp_fasta),
                file(snp_tab),
                file(clusters_tab)

    output:
        path("fasta/*"),        optional: true
        path("positions/*"),    optional: true
        path("nexus/*"),        optional: true

        tuple val(clusterID),
            path("fasta/${clusterID}_refseq.fasta"), emit: snp_fasta

        tuple val(lineage), 
                val(clusterID),
                file("fasta/${clusterID}_refseq.fasta"),
                file(clusters_tab),     emit: annotated_nexus_ch

    script:

        """
        # Get the list of genomes in the cluster
            grep -w "${clusterID}" ${clusters_tab} \\
                | cut -f3 > tmp.sampleIDs

            grep -f tmp.sampleIDs ${snp_fasta} \\
                | sed  's/>//g' \\
                | sort | uniq >genomes.list
        
        seqkit grep -w 0 -f "genomes.list" "${snp_fasta}" > "${clusterID}.fasta"

        # run the nexus script
            bash ${params.script_dir}/shell/create-variable-region-nexus.sh \\
                    -c ${clusterID} \\
                    -f ${clusterID}.fasta \\
                    -t ${snp_tab} \\
                    -m ${params.mtbc_ancestor_path} \\
                    1>>.command.out \\
                    2>>.command.err || true
                    
        # simplify the name of the file
            cat ${snp_tab} > ${clusterID}.snp.tab

        # check if the nexus generation was successful
            if [[ ! -f fasta/${clusterID}_refseq.fasta ]]; then echo "Nexus generation failed" > fasta/${clusterID}_refseq.fasta; fi
            if [[ ! -f positions/${clusterID}_genomic_positions.tab ]]; then echo "Nexus generation failed" > positions/${clusterID}_genomic_positions.tab; fi

        """

}