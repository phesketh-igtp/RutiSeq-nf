process DB_COMPLIANCE_CHECK {
    tag "$params.runID"
    
    publishDir "${params.outDir}/db/", mode: 'copy', pattern: "*.{txt,log}"
    
    input:
        val(sampleID_list)
    
    output:
        val("db_integrity_report.txt"), emit: db_compliance_check
    
    script:
    """
    #!/usr/bin/env Rscript

    #··············································································#
    # Load Environment
    #··············································································#

    # load libraries
        packages <- c( "dplyr", "tidyverse" )

    # Identify missing packages
        missing_pkgs <- packages[!packages %in% installed.packages()[, "Package"]]

    # Install missing packages
        if (length(missing_pkgs) > 0) {
            install.packages(missing_pkgs, dependencies = TRUE)
        }

    # Load all packages
        invisible(lapply(packages, library, character.only = TRUE))

    set.seed(1234)


    #··············································································#
    # Get the database IDs
    #··············································································#
    
    path <- "${params.outDir}/db/samples/"

    dirs <- list.dirs(
                    path, 
                    full.names = FALSE, 
                    recursive = TRUE
                    )

    #··············································································#
    # Check the integrity of the databases
    #··············································································#

    # Get sample IDs (first path element)
    sampleID <- sub("/.*", "", dirs)

    # Create a data frame of unique samples
    df <- data.frame(sampleID = unique(sampleID), stringsAsFactors = FALSE)

    # Helper function to test presence of a subdirectory
    has_dir <- function(sample, dir) {
    any(grepl(paste0("^", sample, "/", dir, "\$"), dirs))
    }

    # Add columns
    df\$mtbseq      <- sapply(df\$sampleID, has_dir, dir = "mtbseq")
    df\$snippy      <- sapply(df\$sampleID, has_dir, dir = "snippy")
    df\$tbprofiler  <- sapply(df\$sampleID, has_dir, dir = "tbprofiler")

    # Finally, check the "name" compliance
    # Count underscores
    underscore_count <- lengths(regmatches(df\$sampleID, gregexpr("_", df\$sampleID)))

    # Error if more than one underscore
    if (any(underscore_count > 1, na.rm = TRUE)) {
            bad <- df\$sampleID[underscore_count > 1]
                stop(
                    paste(
                    "Invalid sampleID(s) with more than one underscore:",
                    paste(bad, collapse = ", ")
                    ),
                    call. = FALSE
            )
        }

    # Ensure there are no empty sampleID fields
    df <- df |> filter(sampleID != "")

    # check all colums are true
    df\$passes <- with(df, mtbseq & snippy & tbprofiler & name_compliant)
    failures <- df |> filter(passes == FALSE)

    # Write report
    write.table(
                df,
                file = "db_integrity_report.txt",
                sep = "\\t",
                row.names = FALSE,
                quote = FALSE
                )

    #··············································································#
    # Halt workflow if any failures
    #··············································································#

    if (nrow(failures) > 0) {
        log_msg <- paste0(
            "DB compliance check failed for ",
            nrow(failures),
            " samples. See db_integrity_report.txt for details."
        )

        write(log_msg, file = "db_compliance_check.log")

        stop(log_msg, call. = FALSE)
    } else {
        write("DB compliance check passed for all samples.", file = "db_compliance_check.log")
    }

    cat("DB compliance check completed successfully.\\n")
    """
}