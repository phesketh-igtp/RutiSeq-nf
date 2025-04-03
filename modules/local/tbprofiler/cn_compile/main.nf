process CN_TBPROFILE_COMPILE {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 0.1
    @description: 
        This module runs TB-Profiler cimpile on the single sample using the TBDB database. 
        It first creates a symbolic link to the data.
*/

    tag "${runID}"
        
    conda params.tbprofiler_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
            else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
            else return null
    }
        
    publishDir "${params.outdir}/bbdd/negative-controls/tbprofiler/", mode: 'copy'

    input:
        val(runID)
        path(all_tbprofiler_results)


    output:
        path("tbprofiler.txt"),         emit: tbprofile_compiled
        path("tbprofiler.variants.csv")
        path("tbprofiler.variants.txt")

    script:
    """
    # Create sybolic links to the tbprofiler results
        ln -s ${params.outdir}/bbdd/negative-controls/tbprofiler/results/ .
        ln -s ${params.outdir}/bbdd/negative-controls/tbprofiler/bam/ .
        ln -s ${params.outdir}/bbdd/negative-controls/tbprofiler/vcf/ .

    # Example command to compile input files
        tb-profiler collate
    """
}