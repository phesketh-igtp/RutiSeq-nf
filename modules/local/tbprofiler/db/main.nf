process TBPROFILER_DB{

    conda "bioconda::tb-profiler=6.3.0"
    
    container 'oras://community.wave.seqera.io/library/tb-profiler:6.3.0--4f362e6be5d39a05'

    output:
    val true

    shell:
    """
    tb-profiler update_tbdb --branch who
    """

}