//
// Subworkflow: Prepare and validate samples from samplesheet
//

include { FETCH_SRA }       from '../modules/prepare-samples_wf/fetch_sra/main'
include { FILE_CHECK }      from '../modules/prepare-samples_wf/init-file-checks/main'

workflow PREPARE_SAMPLES_WF {
    take:
    samplesheet     // path: sample sheet CSV file

    main:

    def color_green  = '\u001B[32m'
    def color_red    = '\u001B[31m'
    def color_reset  = '\u001B[0m'
    def color_cyan   = '\u001B[36m'

    /* 
    // Prepare sample channels
    */

    // Parse sample sheet and validate
    Channel
        .fromPath(samplesheet)
        .ifEmpty { error "Sample sheet file '${samplesheet}' not found or empty" }
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
        .branch {
            sra_ena: it[3] in ['SRA', 'ENA']
            regular: it[3] in ['sample']
            control: it[3] == 'control'
        }
        .set { samples_by_type }

    // Fetch SRA/ENA samples
    FETCH_SRA(
        samples_by_type.sra_ena.map { tuple(it[0], it[3]) }
    )

    // Combine all samples (regular + fetched SRA/ENA)
    samples_ch = samples_by_type.regular
        .map { tuple(it[0], it[1], it[2], it[3]) }  // Include type field
        .mix(
            FETCH_SRA.out.fetch_fastq_tuple
        )

    // Extract controls
    controls_ch = samples_by_type.control
        .map { tuple(it[0], it[1], it[2], it[3]) }  // Include type field here too

    // Check if samples have been previously analyzed
    FILE_CHECK(samples_ch)

    // Collect verified samples
    verified_samples_ch = FILE_CHECK.out.sample_paths
        .collectFile(name: 'all_sample_paths.txt', newLine: true, storeDir: params.outDir)
        .ifEmpty { file("${params.outDir}/empty_all_sample_paths.txt") }

    // Parse verified samples into complete tuple structure
    all_samples_ch = verified_samples_ch
        .splitCsv()
        .map { row -> 
            log.debug "DEBUG - Processing sample row: ${row}, size: ${row.size()}"
            
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
                    mtbseq_vcf && mtbseq_vcf != 'null' && mtbseq_vcf != '' ? file(mtbseq_vcf.trim()) : []
                )   
            } else {
                log.warn "Unexpected row format (size: ${row.size()}): ${row}"
                null
            }
        }
        .filter { it != null }

    // DEBUG: Report final samples
    //all_samples_ch.view { "Final sample: ${it[0]}" }

    emit:
        all_samples   = all_samples_ch     // tuple: [ sampleID, forward, reverse, type, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf ]
        samples       = samples_ch              // tuple: [ sampleID, forward, reverse ]
        controls      = controls_ch            // tuple: [ sampleID, forward, reverse ]

}

/*
    @author: Poppy J Hesketh Best
    @date: 2025-11-12
    @version: 1.0.3
    @description: 
        This workflow prepares and validates samples from a samplesheet. It fetches SRA/ENA samples if needed,
        checks for required columns and non-empty file paths, and verifies if samples have been previously analyzed.
        For complete samples (with existing analysis results), FASTQs remain as empty placeholders.
    @changelog
        v1.0.0-2025-11-12: Initial version
        v1.0.1-2025-11-14: Cleaned up redundancy, FETCH_SRA called once, complete samples keep empty FASTQs
        v1.0.2-2025-11-14: Added 'samples' and 'controls' to emit block
        v1.0.3-2025-11-14: Fixed CSV parsing to handle 11-12 element rows and 'null' strings
*/
