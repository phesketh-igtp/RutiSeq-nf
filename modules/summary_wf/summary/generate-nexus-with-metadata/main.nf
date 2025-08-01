process GENERATE_ANNOTATED_NEXUS {

    conda params.r_stats_env

    tag "cluster: ${clusterID}"

    publishDir "${params.outDir}/results/${params.runID}/networks/nexus/", mode: 'copy'

    input:
        tuple val(clusterID), 
                path(snp_fasta)
        path(metadata)

    output:
        path("${clusterID}.annoated.nex")

    script:
        """

        grep -f genomes.list ${metadata} | cut -d ',' -f2,3 > dates.csv

        bash ${params.script_dir}/shell/add-nexus-dates.sh \\
            -i ${snp_fasta} -m dates.csv -p ${clusterID}

        """

}