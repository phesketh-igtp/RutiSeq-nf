//
// Subworkflow: Prepare and validate samples from samplesheet
//
include { FETCH_SRA }       from '../modules/prepare-samples_wf/fetch_sra/main'
include { FILE_CHECK }      from '../modules/prepare-samples_wf/init-file-checks/main'
//include { VERSION_LOGGING } from '../modules/version-logging/main'

workflow PREPARE_SAMPLES_WF {
    take:
    samplesheet     // path: sample sheet CSV file

    main:

    def green  = '\u001B[32m'
    def red    = '\u001B[31m'
    def reset  = '\u001B[0m'
    def cyan   = '\u001B[36m'

    /* 
        Prepare sample channels
    */

    // Parse sample sheet and validate samplesheet file exists
    def samplesheet_file = file(samplesheet)
    if (!samplesheet_file.exists()) {
        error "${red}Sample sheet file ${cyan}'${samplesheet}'${red} not found${reset}"
    }
    if (samplesheet_file.isEmpty()) {
        error "${red}Sample sheet file ${green}'${samplesheet}'${red} is empty${reset}"
    }

    // Parse sample sheet and validate
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def requiredColumns = ['originalID', 'sampleID', 'forward_path', 'reverse_path', 'type']
            def missingColumns = requiredColumns.findAll { !row.containsKey(it) }
            if (missingColumns) {
                error "Missing required column(s) in samplesheet: ${missingColumns.join(', ')}"
            }
            
            def sampleType = row.type.trim()
            
            // Validate paths for non-SRA/ENA samples
            if (sampleType !in ['SRA', 'ENA']) {
                if (!row.forward_path.trim() || !row.reverse_path.trim()) {
                    error "Empty file path found for sample ${row.sampleID}. Both forward and reverse paths must be provided."
                }
            }
            
            // Create file objects or empty lists for SRA/ENA
            def forwardFile = (sampleType in ['SRA', 'ENA']) ? [] : file(row.forward_path.trim(), checkIfExists: false)
            def reverseFile = (sampleType in ['SRA', 'ENA']) ? [] : file(row.reverse_path.trim(), checkIfExists: false)
            
            tuple(row.sampleID.trim(), forwardFile, reverseFile, sampleType)
        }
        .set { parsed_samples_ch }

    // STEP 1: First check all samples (including SRA/ENA) for existing analysis files
    FILE_CHECK(parsed_samples_ch)

    // STEP 2: Parse FILE_CHECK output and determine which samples need SRA download
    FILE_CHECK.out.sample_paths
        .map { row -> 
            log.debug "DEBUG - Processing FILE_CHECK output row: ${row}, size: ${row.size()}"
            
            // Handle rows with 11 or 12 elements (with trailing commas/empty fields)
            if (row.size() >= 11) {
                def sampleID = row[0]
                def forward = row[1]
                def reverse = row[2]
                def type = row[3]
                def mtbseq_class = row[4]
                def mtbseq_stats = row[5]
                def mtbseq_pos = row[6]
                def mtbseq_vars = row[7]
                def tbdb_out = row[8]
                def who_out = row[9]
                def mtbseq_vcf = row[10]
                
                // Check if all analysis files are empty (indicating need for SRA download)
                def analysisFilesEmpty = [mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf]
                    .every { it == null || it == '' || it == 'null' }
                
                // Check if FASTQ files are also empty (for SRA/ENA samples)
                def fastqFilesEmpty = (forward == null || forward == '' || forward == 'null') && 
                                        (reverse == null || reverse == '' || reverse == 'null')
                
                tuple(
                    sampleID,
                    forward && forward != 'null' && forward != '' ? file(forward.trim()) : [],
                    reverse && reverse != 'null' && reverse != '' ? file(reverse.trim()) : [],
                    type && type != 'null' && type != '' ? type : 'sample',
                    mtbseq_class && mtbseq_class != 'null' && mtbseq_class != '' ? file(mtbseq_class.trim()) : [],
                    mtbseq_stats && mtbseq_stats != 'null' && mtbseq_stats != '' ? file(mtbseq_stats.trim()) : [],
                    mtbseq_pos && mtbseq_pos != 'null' && mtbseq_pos != '' ? file(mtbseq_pos.trim()) : [],
                    mtbseq_vars && mtbseq_vars != 'null' && mtbseq_vars != '' ? file(mtbseq_vars.trim()) : [],
                    tbdb_out && tbdb_out != 'null' && tbdb_out != '' ? file(tbdb_out.trim()) : [],
                    who_out && who_out != 'null' && who_out != '' ? file(who_out.trim()) : [],
                    mtbseq_vcf && mtbseq_vcf != 'null' && mtbseq_vcf != '' ? file(mtbseq_vcf.trim()) : [],
                    analysisFilesEmpty && fastqFilesEmpty && (type in ['SRA', 'ENA']) // needs_sra_download flag
                )   
            } else {
                log.warn "Unexpected row format (size: ${row.size()}): ${row}"
                null
            }
        }
        .filter { it != null }
        .branch {
            needs_sra: it[11] == true  // Last element is the needs_sra_download flag
            complete: it[11] == false  // Sample has existing files or is not SRA/ENA
        }
        .set { checked_samples }

    // STEP 3: Prepare SRA samples for FETCH_SRA (with full tuple structure)
    sra_samples_to_fetch = checked_samples.needs_sra
        .map { sample ->
            // Create the full tuple that FETCH_SRA expects
            tuple(
                sample[0],  // accession (sampleID)
                sample[1],  // forward (empty for SRA)
                sample[2],  // reverse (empty for SRA)
                sample[3],  // type
                sample[4],  // mtbseq_class
                sample[5],  // mtbseq_stats
                sample[6],  // mtbseq_pos
                sample[7],  // mtbseq_vars
                sample[8],  // tbdb_out
                sample[9],  // who_out
                sample[10]  // mtbseq_vcf
            )
        }

    // STEP 4: Download SRA files only for samples that need them
    FETCH_SRA(sra_samples_to_fetch)

    // STEP 5: Combine complete samples with fetched SRA samples
    all_samples_ch = checked_samples.complete
        .map { tuple(it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7], it[8], it[9], it[10]) }  // Remove the flag
        .mix(FETCH_SRA.out.fetch_fastq_tuple)

    // Create simplified channels for backward compatibility
    samples_ch = all_samples_ch
        .filter { it[3] == 'sample' }
        .map { tuple(it[0], it[1], it[2]) }  // sampleID, forward, reverse

    controls_ch = all_samples_ch
        .filter { it[3] == 'control' }
        .map { tuple(it[0], it[1], it[2]) }  // sampleID, forward, reverse

    // Collect all samples for output file
    verified_samples_ch = all_samples_ch
        .map { sample ->
            // Convert tuple back to CSV format for collectFile
            "${sample[0]},${sample[1] ?: ''},${sample[2] ?: ''},${sample[3]},${sample[4] ?: ''},${sample[5] ?: ''},${sample[6] ?: ''},${sample[7] ?: ''},${sample[8] ?: ''},${sample[9] ?: ''},${sample[10] ?: ''}"
        }
        .collectFile(name: 'all_sample_paths.txt', newLine: true, storeDir: params.outDir)
        .ifEmpty { file("${params.outDir}/empty_all_sample_paths.txt") }

    emit:
        all_samples   = all_samples_ch     // tuple: [ sampleID, forward, reverse, type, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf ]
        samples       = samples_ch         // tuple: [ sampleID, forward, reverse ]
        controls      = controls_ch        // tuple: [ sampleID, forward, reverse ]
}

/*
@author: Poppy J Hesketh Best
@date: 2025-11-12
@version: 1.0.7
@description: 
    This workflow prepares and validates samples from a samplesheet. It first checks all samples for existing analysis files.
    Only samples with empty analysis indexes (excluding sampleID) and that are SRA/ENA type will have their FASTQ files downloaded.
    The FETCH_SRA module receives the full tuple structure and returns updated tuples with downloaded FASTQ paths.
@changelog
    v1.0.0-2025-11-12: Initial version
    v1.0.1-2025-11-14: Cleaned up redundancy, FETCH_SRA called once, complete samples keep empty FASTQs
    v1.0.2-2025-11-14: Added 'samples' and 'controls' to emit block
    v1.0.3-2025-11-14: Fixed CSV parsing to handle 11-12 element rows and 'null' strings
    v1.0.4-2026-01-05: Improved error handling with samplesheet
    v1.0.5-2026-01-19: Modified to first check files, then download SRA samples for better workflow efficiency
    v1.0.6-2026-01-19: Conditional SRA download only for samples with empty analysis indexes
    v1.0.7-2026-01-19: Fixed FETCH_SRA input to match module's expected tuple structure
*/