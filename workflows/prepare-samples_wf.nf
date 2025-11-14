//
// Subworkflow: Prepare and validate samples from samplesheet
//

include { FETCH_SRA }   from '../modules/prepare-samples_wf/fetch_sra/main'
include { FILE_CHECK }  from '../modules/prepare-samples_wf/init-file-checks/main'

workflow PREPARE_SAMPLES {
    take:
    samplesheet     // path: sample sheet CSV file

    main:

    def color_green  = '\u001B[32m'
    def color_red    = '\u001B[31m'
    def color_reset  = '\u001B[0m'
    def color_cyan   = '\u001B[36m'

    // Create channel from sample sheet

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
            
            // Only check for empty paths if NOT SRA/ENA (they will be downloaded later)
            if (sampleType !in ['SRA', 'ENA']) {
                if (!row.forward_path.trim() || !row.reverse_path.trim()) {
                    error "Empty file path found for sample ${row.sampleID}. Both forward and reverse paths must be provided."
                }
            }
            
            // For SRA/ENA, create placeholder empty lists; for others, use file paths
            def forwardFile = (sampleType in ['SRA', 'ENA']) ? [] : file(row.forward_path.trim(), checkIfExists: false)
            def reverseFile = (sampleType in ['SRA', 'ENA']) ? [] : file(row.reverse_path.trim(), checkIfExists: false)
            
            tuple(
                row.sampleID.trim(), 
                forwardFile, 
                reverseFile, 
                sampleType
            )
        }
        .branch {
            sample: it[3] in ['sample', 'SRA', 'ENA']
            control: it[3] == 'control'
        }
        .set { branched_samples_by_type }

    // Further split samples into SRA/ENA vs regular samples
    branched_samples_by_type.sample
        .branch {
            sra_ena: it[3] in ['SRA', 'ENA']
            regular: true
        }
        .set { sample_types }

    // Fetch SRA/ENA samples
    FETCH_SRA(
        sample_types.sra_ena.map { tuple(it[0], it[3]) }
    )

    // Update SRA/ENA samples with downloaded FASTQs
    updated_sra_samples = FETCH_SRA.out.fetch_fastq_tuple
        .map { sampleID, forward, reverse ->
            tuple(sampleID, forward, reverse)
        }

    // Merge regular samples with updated SRA/ENA samples
    samples_ch = sample_types.regular
        .map { tuple(it[0], it[1], it[2]) }
        .mix(updated_sra_samples)

    // Process controls (remove type field)
    controls_ch = branched_samples_by_type.control
        .map { tuple(it[0], it[1], it[2]) }

    // Report the samples
    samples_ch.view { sampleID, forward, reverse ->
        "${color_cyan}Sample: ${color_green}${sampleID}${color_reset} | ${color_cyan}Forward: ${color_green}${forward}${color_reset} | ${color_cyan}Reverse: ${color_green}${reverse}${color_reset}"
    }

    // Report the controls
    controls_ch.view { sampleID, forward, reverse ->
        "${color_red}Control: ${color_green}${sampleID}${color_reset} | ${color_red}Forward: ${color_green}${forward}${color_reset} | ${color_red}Reverse: ${color_green}${reverse}${color_reset}"
    }

    // Check if genomes have previously been analyzed
    FILE_CHECK(samples_ch)

    // Process verified samples
    verified_samples_ch = FILE_CHECK.out.sample_paths
        .collectFile(name: 'all_sample_paths.txt', newLine: true, storeDir: params.outDir)
        .ifEmpty { file("${params.outDir}/empty_all_sample_paths.txt") }

    // Parse samples into desired tuple structure
    comp_samples_ch = verified_samples_ch
        .splitCsv()
        .map { row -> 
            log.debug "DEBUG - Processing sample row: ${row}"
            if (row.size() == 10) {
                def (sampleID, forward, reverse, type, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = row
                tuple(
                    sampleID,
                    forward ? file(forward.trim()) : [],
                    reverse ? file(reverse.trim()) : [],
                    type,
                    mtbseq_class ? file(mtbseq_class.trim()) : [],
                    mtbseq_stats ? file(mtbseq_stats.trim()) : [],
                    mtbseq_pos ? file(mtbseq_pos.trim()) : [],
                    mtbseq_vars ? file(mtbseq_vars.trim()) : [],
                    tbdb_out ? file(tbdb_out.trim()) : [],
                    who_out ? file(who_out.trim()) : [],
                    mtbseq_vcf ? file(mtbseq_vcf.trim()) : []
                )   
            } else {
                log.warn "Error with channel: ${row}"
                null
            }
        }
        .filter { it != null }

    // Branch the complete samples for SRA/ENA processing
    comp_samples_ch
        .branch {
            sra_ena: it[3] in ['SRA', 'ENA']
            other: true
        }
        .set { branched_ENA_samples }

    // Fetch SRA/ENA samples (if any in comp_samples_ch)
    FETCH_SRA(
                branched_ENA_samples.sra_ena.map { tuple(it[0], it[3]) }
            )

    // Update the tuple with downloaded FASTQs
    updated_comp_sra_samples = branched_ENA_samples.sra_ena
        .map { tuple ->
            def sampleID = tuple[0]
            def type = tuple[3]
            def mtbseq_class = tuple[4]
            def mtbseq_stats = tuple[5]
            def mtbseq_pos = tuple[6]
            def mtbseq_vars = tuple[7]
            def tbdb_out = tuple[8]
            def who_out = tuple[9]
            def mtbseq_vcf = tuple[10]
            
            tuple(sampleID, type, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf)
        }
        .join(FETCH_SRA.out.fetch_fastq_tuple)
        .map { tuple ->
            def sampleID = tuple[0]
            def type = tuple[1]
            def mtbseq_class = tuple[2]
            def mtbseq_stats = tuple[3]
            def mtbseq_pos = tuple[4]
            def mtbseq_vars = tuple[5]
            def tbdb_out = tuple[6]
            def who_out = tuple[7]
            def mtbseq_vcf = tuple[8]
            def forward_fastq = tuple[9]
            def reverse_fastq = tuple[10]
            
            tuple(sampleID, forward_fastq, reverse_fastq, type, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf)
        }

    // Merge back with non-SRA/ENA samples
    all_samples_ch = branched_ENA_samples.other.mix(updated_comp_sra_samples)

    // Final channel view
    all_samples_ch.view { "Final sample: ${it[0]}" }

    emit:
        all_samples = all_samples_ch     // tuple: [ sampleID, forward, reverse, type, ... all other fields ]
    
}

/*
    @author: Poppy J Hesketh Best
    @date: 2025-11-12
    @version: 1.0.0
    @description: 
        This workflow prepares and validates samples from a samplesheet. It fetches SRA/ENA samples if needed,
        checks for required columns and non-empty file paths, and verifies if samples have been previously analyzed.
        Was originally in the main.nf, but moved to its own workflow for better modularity.
    @changelog
        v1.0.0-2025-11-12: Initial version
*/