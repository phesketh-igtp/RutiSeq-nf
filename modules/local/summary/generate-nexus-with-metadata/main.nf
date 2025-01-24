process GENERATE_NEXUS_W_METADATA {

    conda params.snp_profiling_env 

    tag "cluster: ${clusterID}"

    publishDir "${params.outdir}/results/networks/nexus/", mode: 'copy'

    input:
        tuple val(clusterID), 
                path(nexus)
        path(metadata)
        path(pairwise_clusters)

    output:
        path("nexus/${clusterID}.time.nex")
        path("nexus/${clusterID}.loc.nex")

    script:
        """
        pairwise_clusters |> genomes.list

        grep -f genomes.list ${metadata} | cut -f2,3 > dates.csv
        grep -f genomes.list ${metadata} | cut -f2,4 > loc.csv

        Rscript -e "library(tidyverse)
        dates <- read.delim('dates.csv', header = TRUE, sep = ",") |>
                    mutate(seen = 1 ) |>
                    pivot_wider(names_from = location, values_from = seen)
        dates[is.na(dates)] <- 0

        loc <- read.delim('loc.csv', header = TRUE, sep = ",") |>
                    mutate(seen = 1 ) |>
                    pivot_wider(names_from = location, values_from = seen)
        loc[is.na(loc)] <- 0

        write.table(dates)
        write.table(loc)

        "

        metadata |> location 

        """

}