include { INSPECT_BBDD               }  from '../modules/local/negative-ctrls/inspect_bbdd/main.nf'
include { CN_READ_TAXONOMY           }  from '../modules/local/negative-ctrls/inspect_reads/main.nf'
include { CN_TBPROFILER_TBDB         }  from '../modules/local/tbprofiler/cn_profile.tbdb/main.nf'
include { CN_MTBSEQ_SINGLE           }  from '../modules/local/mtbseq/cn_single/main.nf'
//include { COMPILE_CN_READS_SUMMARY   }  from '../modules/local/negative-ctrls/combine-qc-results/main.nf'

workflow NEGATIVE_CTRL_WF {

    take:
        controls_ch
        tbprofiler_update_db

    main:

        /*
        Run KAIJU on the reads and get read taxonomy
        */
            INSPECT_BBDD(controls_ch)

            // After the FILE_CHECK process
            verified_controls_ch = INSPECT_BBDD.out.controls_paths
                .collectFile(name: 'all_controls_paths.txt', newLine: true, storeDir: params.outdir)
                .ifEmpty { file("${params.outdir}/empty_all_controls_paths.txt") }

            // Parse the controls into the desired tuple structure
                comp_controls_ch = verified_controls_ch
                    .splitCsv()
                    .map { row -> 
                        log.debug "DEBUG - Processing sample row: $row"
                        if (row.size() == 3) {
                            def (sampleID, forward, reverse) = row
                            tuple(
                                sampleID,
                                forward ? file(forward.trim()) : [],
                                reverse ? file(reverse.trim()) : [],
                            )   
                        } else {
                            log.warn "Error with channel: $row"
                            null
                        }
                    }
                    .filter { it != null }

                // Demonstrate the content of the channel
                /// comp_controls_ch.view { sample -> "Sample: $sample" }

            // Branch the channel into those with outputs and those without:
                branched_channel = comp_controls_ch.branch {
                with_reads: it[1] != [] && it[2] != [] // zero-indexed so [1] is the second value in the tuple, ect
                without_reads: it[1] == [] || it[2] == [] }

        control_ch_analysis = branched_channel.with_reads

        /*
        Run KAIJU on the reads and get read taxonomy
        */

            CN_READ_TAXONOMY( control_ch_analysis )

        // collect all the results
            //all_cn_k2_results   = CN_READ_TAXONOMY.out.cn_k2_results.map { it[1] }.collect()
            //all_cn_stats        = CN_READ_TAXONOMY.out.cn_stats.map { it[1] }.collect()

            //all_cn_k2_results.view()
            //all_cn_stats.view()

        /*
            Run Tb-Profiler and MTBseq on the reads (expect them to fail)
        */
            CN_TBPROFILER_TBDB( control_ch_analysis, tbprofiler_update_db )
            CN_MTBSEQ_SINGLE( control_ch_analysis )

        /*
            Compile the Negative control read summary
        */
            //COMPILE_CN_READS_SUMMARY(all_cn_k2_results, all_cn_stats)

}