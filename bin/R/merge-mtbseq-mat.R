# Load necessary library
library(optparse)

# Define command-line options
    option_list <- list(
        make_option(
        c("-d", "--directory"), 
        type = "character", 
        default = NULL, 
        help = "Path to the directory containing TSV files", 
        metavar = "character"
        )
    )

# Parse the command-line options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

# Check if directory is provided
    if (is.null(opt$directory)) {
        stop("Please provide a directory path using the -d or --directory flag.")
    }

# Get the directory path
    directory_path <- opt$directory

# Verify that the directory exists
    if (!dir.exists(directory_path)) {
        stop("The provided directory does not exist.")
    }

# List all TSV files in the directory
    file_list <- list.files(path = directory_path, pattern = "\\.tsv$", full.names = TRUE)

# Extract the base names (without extensions) and identify duplicates
    base_names <- gsub("_\\d+\\.tsv$", "", basename(file_list))
    unique_files <- file_list[!duplicated(base_names)]

# Resolve duplicates by keeping the file with the highest numeric value
    highest_files <- sapply(unique(base_names), function(name) {
        matching_files <- grep(paste0("^", name, "_\\\\d+\\.tsv$"), basename(file_list), value = TRUE)
        numeric_values <- as.numeric(gsub(".*_(\\\\d+)\\.tsv$", "\\1", matching_files))[which.max(numeric_values)]
        }
    )

# Full paths to the highest value files
    file_list <- file.path(directory_path, highest_files)

# Check if there are any TSV files
    if (length(file_list) == 0) {
        stop("No TSV files found in the specified directory.")
    }

# Import all TSV files into a list of data frames
    data_list <- lapply(file_list, function(file) {
        read.delim(file, header = TRUE, row.names = 1)
        }
    )

# Optionally, name the list elements based on file names
    names(data_list) <- basename(file_list)

# Merge all data frames into a single matrix
    merged_matrix <- do.call(cbind, data_list)

# Print the summary of imported and merged data
    cat("Successfully imported and merged the following TSV files into a single matrix:\n")
    cat(paste(basename(file_list), collapse = "\n"), "\n")

# Optionally: Save the merged matrix to a file
    write.table(merged_matrix, file = "Master-matrixes.tsv", sep = "\t", quote = FALSE)
