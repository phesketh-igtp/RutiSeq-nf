process TBPROFILER_DB{

    conda "bioconda::tb-profiler=6.3.0"
    
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'

    output:
    val true

    shell:
    """
    tb-profiler update_tbdb --branch who
    """

}