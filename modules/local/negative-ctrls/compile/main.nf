process CN_READS_SUMMARY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-03
    @version: 0.1
    @description:
        This module concatenates all the statistical files within the output directory
        capturing all the results to date and not just in this run
                
*/

    tag "${runID}"
        
    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
            else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
            else return null
    }
        
    publishDir "${params.outdir}/bbdd/negative-controls/mtbseq/", mode: 'copy'

    input:
        val(runID)
        path(all_cn_k2_results)
        path(all_cn_stats)


    output:
        path("negative-controls.k2.report"),    emit: k2_combined
        path("negative-controls.stats.tsv"),    emit: stats_combined

    script:
    """
    # Concatenate the k2.results
        for file in ${params.outdir}/bbdd/negative-controls/Classification/*.k2.report; do
            cat \${file} | sed '1d' >> negative-controls.k2.report
        done

    # Concatenate the read statistics
        for file in ${params.outdir}/bbdd/negative-controls/Statistics/*.stats.tsv; do
            cat \${file} | sed '1d' >> negative-controls.stats.tsv
        done
    """
}