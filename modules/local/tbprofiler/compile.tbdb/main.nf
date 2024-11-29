process TBPROFILER_COMPILE_TBDB {
    
    tag "$params.runID"

    conda (params.enable_conda ? { file("/imppc/labs/emlab/phesketh/miniconda3/envs/tb-profiler").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/tb-profiler" : "./modules/local/tbprofiler/conda.yml" } : null)

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'

    publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'copy'

    input:
    path tbprofiler_tbdb_txt
    path tbprofiler_tbdb_json

    output:
    path "tbprofile.txt",  emit: tbprofile_tbdb_compile
    path "tbprofile..variants.csv" 
    path "tbprofiler.variants.txt"

    script:
    """
    tb-profiler collate --full --mark_missing --all_variants --itol
    """
}