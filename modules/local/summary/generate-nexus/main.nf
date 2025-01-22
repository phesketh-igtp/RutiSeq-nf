process GENERATE_NEXUS {

    conda 

    tag "${lineage}"

    publishDir "${params.outdir}/bbdd/results/networks/", mode: 'copy'

    input:
        path pairwise_clusters
        tuple val(lineage), 
                path(joint_dir), 
                path(amend_dir)

    output:
        path("mj-networks/*")
        path("mj-networks/nexus/*"),            emit: nexus_dir
        path("mj-networks/fasta/*")
        path("mj-networks/positions/*")

    script:
        """
        
        """

}