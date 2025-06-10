process GENERATE_NEXUS_W_METADATA {

    conda params.r_stats_env

    tag "cluster: ${clusterID}"

    publishDir "${params.outdir}/results/networks/nexus/", mode: 'copy'

    input:
        tuple val(clusterID), 
                path(nexus)
        path(metadata)

    output:
        path("${clusterID}_dates.nex")
        path("${clusterID}_locs.nex")

    script:
        """

        grep -f genomes.list ${metadata} | cut -d ',' -f2,3 > dates.csv
        grep -f genomes.list ${metadata} | cut -d ',' -f2,4 > loc.csv

        bash ${params.script_dir}/shell/add-nexus-metadata_dates.sh \\
            -i ${nexus} -p ${clusterID}

        """

}