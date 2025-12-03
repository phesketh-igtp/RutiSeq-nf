process PREPROCESS_CLUSTER {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0.1
    @description:
        This process runs the MTBseq TBgroups step on the joint and amend directories
        for each lineage and distance. It takes the output from the
        MTBSEQ_LINEAGE_JOINT_AMEND() process and runs the TBgroups step on the joint
        and amend directories. It also renames the output files for simplicity.
        It also wrangles the output matrix into a useful format for downstream analysis.
        Occasionally, the MTBseq TBgroups step will fail and produice empty files/no-files.
    @last_updated: 2025-04-01
    @changelog:
        v1.0.0-2025-04-01: Initial version + documnetadocumentation
        v1.0.1-2025-05-19: Removed additonal_args as it was not utilised and was creating inconsistencies
*/

    tag "${lineage}; t=${distance}"

    conda params.r_stats_env

    // TODO: container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
    ///        } else { 'quay.io/biocontainers/mtbseq' }
    ///}
                
    publishDir "${params.outDir}/db/comparison/mtbseq/${lineage}/Groups/", mode: 'copy', pattern: ".clusters.tsv"

    input:
        // Nexus output
        tuple val(lineage), 
            val(distance), 
            path(fasta),
            path(tab),
            path(clusters, stageAs: 'clusters_input.tsv'),
            path(singletons)

    output:
        // Nexus output
        tuple val(lineage), 
            val(distance), 
            path("${lineage}_snps.fasta"),
            path("${lineage}_snps.tab"),
            path("${lineage}_d${distance}.clusters.tsv"), emit: nexus_ch

        path("${lineage}_d${distance}.clusters.tsv"), emit: pairwise_clusters_processed

    script:

        //def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        #!/usr/bin/env Rscript --vanilla

        # Enable comprehensive error reporting
        options(error = function() {
            cat("ERROR occurred at:", date(), "\\n")
            cat("Traceback:\\n")
            traceback()
            quit(status = 1)
        })

        cat("=== CONCATENATE_CLUSTERS DEBUG START ===\\n")
        cat("R version:", R.version.string, "\\n")
        cat("Working directory:", getwd(), "\\n")
        cat("Files in working directory:\\n")
        print(list.files(recursive = TRUE))

        #--------------------------------------------------------------------------#
        # Load packages
        #--------------------------------------------------------------------------#

        # Set a default CRAN mirror (avoids prompt)
        options(repos = c(CRAN = "https://cloud.r-project.org"))

        # List of packages
        packages <- c("tidyverse", "dplyr")

        # Check installed packages
        cat("Installed packages:\\n")
        installed_pkgs <- installed.packages()[, "Package"]
        cat("tidyverse installed:", "tidyverse" %in% installed_pkgs, "\\n")
        cat("dplyr installed:", "dplyr" %in% installed_pkgs, "\\n")

        # Identify missing packages
        missing_pkgs <- packages[!packages %in% installed_pkgs]
        cat("Missing packages:", paste(missing_pkgs, collapse = ", "), "\\n")

        # Install missing packages (non-interactive)
        if (length(missing_pkgs) > 0) {
            cat("Installing missing packages...\\n")
            tryCatch({
                install.packages(missing_pkgs, dependencies = TRUE)
                cat("Package installation completed\\n")
            }, error = function(e) {
                cat("ERROR installing packages:", e\$message, "\\n")
                quit(status = 1)
            })
        }

        # Load all packages
        cat("Loading packages...\\n")
        for (pkg in packages) {
            tryCatch({
                library(pkg, character.only = TRUE)
                cat("Successfully loaded:", pkg, "\\n")
            }, error = function(e) {
                cat("ERROR loading package", pkg, ":", e\$message, "\\n")
                quit(status = 1)
            })
        }

        cat("All packages loaded successfully\\n")

        #--------------------------------------------------------------------------#
        # Read input files
        #--------------------------------------------------------------------------#

        cat("Reading input files...\\n")

        clusters <- read.table(
            file = "clusters_input.tsv",
            sep = "\\t",
            header = FALSE,
            stringsAsFactors = FALSE
        )

        #--------------------------------------------------------------------------#
        # Process clusters
        #--------------------------------------------------------------------------#

        cat("Processing clusters...\\n")

        clusters <- clusters |>
            distinct() |>
            tidyr::separate_wider_delim(
                V4,
                delim = "-",
                names = c("num", "dist", "lin"),
                too_few = "align_start"
            ) |>
            dplyr::mutate(
                num  = stringr::str_pad(num,  width = 3, pad = "0"),
                dist = stringr::str_pad(dist, width = 2, pad = "0")
            ) |> 
            mutate(V4 = paste0(num, "-", 
                                dist, "-", 
                                lin,
                                sep = "")) |>
            select(V1, V2, V3, V4)

        #--------------------------------------------------------------------------#
        # Output file
        #--------------------------------------------------------------------------#
        
        cat("Writing clusters file...\\n")

        write.table(
            clusters,
            file = "${lineage}_d${distance}.clusters.tsv",
            sep = "\\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )
        """
}