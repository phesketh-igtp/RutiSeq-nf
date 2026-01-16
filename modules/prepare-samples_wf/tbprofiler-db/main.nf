process TBPROFILER_DB {

/*
    @author: Poppy J hesketh Best
    @date: 2025-11-17
    @version: v1.0.0
    @description:
        Download TB-Profiler database and sotres in the storeDir, meaning
            if it exists in the storeDir, then it will not download the files.
        Outputs a channel of the database that are used by the respective tools.
    @changelog:
        v1.0.0-2025-11-17: Functioning module created.
*/

    tag "TB-Profiler database: TBDB and WHO"
    
    conda params.tbprofiler_env
    
    storeDir "${params.storeDir}/tbprofiler/"

    input:
        val(runID)

    output:
        path "tbdb/", emit: db
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    # Create the database directory
        mkdir -p tbdb
        
    # Download/update the TB-Profiler database
        tb-profiler update_tbdb
        tb-profiler update_tbdb --branch who
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    runID: ${runID}
        tbprofiler: \$(tb-profiler version 2>&1 | grep "TBProfiler version" | sed 's/.*version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p tbdb
    touch tbdb/tbdb.fasta
    touch tbdb/tbdb.bed
    touch tbdb/tbdb.gff
    touch tbdb/tbdb.ann
    touch tbdb/tbdb.snpeff.config
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    runID: stub
        tbprofiler: \$(tb-profiler version 2>&1 | grep "TBProfiler version" | sed 's/.*version //g')
    END_VERSIONS
    """
}