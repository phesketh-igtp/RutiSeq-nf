process GENERATE_NEXUS_W_MRCA {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process generates a nexus file with the MRCA of the cluster and the H37Rv reference genome.
        It also generates a tab file with the genomic positions of the variants.
        The process takes the output of the SNP profiling process and the pairwise clusters file.
*/

    conda params.snp_profiling_env 

    tag "cluster: ${clusterID}"

    publishDir "${params.outDir}/results/${runID}/networks/", mode: 'copy'

    input:
        val(runID)
        tuple val(lineage), 
                val(clusterID),
                file(snp_fasta),
                file(snp_tab),
                file(ancestor)
                
        file(pairwise_clusters)

    output:
        tuple val(clusterID),
            path("nexus/${clusterID}_refseq_mrca.nex", optional: true), emit: nexus_w_no_metadata
        path("fasta/*"),        optional: true
        path("positions/*"),    optional: true
        path("nexus/*"),        optional: true

    script:

        """
        processed_clusters
        bash ${params.script_dir}/shell/create-variable.region.nexus.w.MRCA.sh \
                -c ${clusterID} \\
                -p ${pairwise_clusters} \\
                -f ${snp_fasta} \\
                -t ${snp_tab} \\ 
                -m ${params.mtbc_ancestor_path} \\
                -n ${ancestor} \\
                -l ${lineage}  \\
                1>>.command.out \\
                2>>.command.err || true

        bash ${params.script_dir}/shell/add-nexus-dates.sh \\
                    -i fasta/${clusterID}_refseq.fasta \\
                    -m ${params.metadata} \\
                    -p nexus/${clusterID}

        """

}
