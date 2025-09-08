process GENERATE_ANNOTATED_NEXUS {

    conda params.r_stats_env

    tag "cluster: ${clusterID}"

    publishDir "${params.outDir}/results/${params.runID}/networks/nexus/", mode: 'copy'

    input:
        tuple val(lineage), 
                val(clusterID),
                file(snp_fasta),
                file(clusters_tab)
        path(metadata)

    output:
        path("${clusterID}.annoated.nex", optional: true)

    script:
        """
        # Get the list of genomes in the cluster
            grep -w "${clusterID}" ${clusters_tab} \\
                | cut -f3 \\
                | grep - ${snp_fasta} \\
                | sed  's/>//g' \\
                | sort | uniq > genomes.list

        if [ ! -s genomes.list ]; then

            echo "The genomes list is empty. Not creating "
        
        else

            # Filter the metadata for just the dates and sampleIDs
            grep -f genomes.list ${metadata} | cut -d ',' -f2,3 > dates.csv

            bash ${params.script_dir}/shell/add-nexus-dates.sh \\
                -i ${snp_fasta} -m dates.csv -p ${clusterID}

        fi  
        """

}