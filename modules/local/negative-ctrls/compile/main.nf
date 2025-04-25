process CN_READS_SUMMARY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-03
    @version: 1.2.0
    @description:
        This module concatenates all the statistical files within the output directory
        capturing all the results to date and not just in this run
    @changelog:
        v1.0.0-2025-04-03:
            Added - initial version
        v1.1.0-2025-04-04:
            Added - taxonkit to get formatted taxonomy
            Changed - output from TSV to CSV
        v1.2.0-2025-04-10:
            Changed - from 'taxonkit lineage' to 'taxonkit reformat' for consistency
*/

    tag "${runID}"

    conda params.taxonkit_env
        
    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_taxonkit
            else if (workflow.containerEngine == 'docker') return params.docker_taxonkit
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_taxonkit
    }
        
    publishDir "${params.outDir}/negative-controls/", mode: 'copy'

    input:
        val(runID)
        val(complete_sampleID)
        path(taxonkit_update_db)


    output:
        path("negative-controls.k2.report")
        path("negative-controls.k2.report.csv"), emit: k2_combined
        path("negative-controls.stats.csv"),    emit: stats_combined

    script:
    """
    # Concatenate the k2.results
        for file in ${params.outDir}/negative-controls/Classification/*.k2.report; do
            cat \${file} | sed '1d' | cut -f1,2,3,6 >> negative-controls.k2.report
        done

    # Run taxonkit to get nicely formatted taxonomy
        taxonkit reformat \\
            -f "{p};{c};{o};{f};{g};{s};{T}" \\
            -F -t -i 4 \\
            --data-dir ${params.outDir}/db/taxonkit/ \\
            negative-controls.k2.report \\
            > negative-controls.k2.report.tmp

    # Create header for final output
        echo "sampleID,percentage,num_reads,taxID,taxonomy" > negative-controls.k2.report.csv

    # Wrangle the output into the correct format
        sed 's@\t@,@g' negative-controls.k2.report.tmp > negative-controls.k2.report.tmp2
        cat negative-controls.k2.report.tmp2 >> negative-controls.k2.report.csv
        rm negative-controls.k2.report.tmp2 negative-controls.k2.report.tmp

    # Concatenate the read statistics
        echo "sampleID,file,format,type,num_seqs,sum_le,min_len,avg_len,max_len,Q1,Q2,Q3,sum_gap,N50,N50_num,Q20(%),Q30(%),AvgQual,GC(%),sum_n" > negative-controls.stats.csv

        for file in ${params.outDir}/negative-controls/Statistics/*.stats.tsv; do
            cat \${file} | sed '1d' | sed 's@\t@,@g' >> negative-controls.stats.csv
        done
    """
}