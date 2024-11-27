process CHECK_EXISTING_OUTPUTS {

    tag "${sampleID}"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    
    output:
    tuple val(sampleID), path(forward), path(reverse), val(tbprofiler_tbdb_exists), val(tbprofiler_who_exists), val(mtbseq_exists)

    script:
    
    """
    tbprofiler_tbdb_exists=0
    tbprofiler_who_exists=0
    mtbseq_exists=0

    if [ -f "${params.outdir}/bbdd/tbprofiler/${sampleID}.results.json" && "${params.outdir}/TBPROFILER_PROFILE_TBDB/${sampleID}.txt" ]; then
        tbprofiler_tbdb_exists=1
    fi

    if [ -f "${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}.results.json" ]; then
        tbprofiler_who_exists=1
    fi

    if [ -f "${params.outdir}/bbdd/mtbseq/${sampleID}/Classification/Strain_Classification.tab" && "${params.outdir}/bbdd/mtbseq/${sampleID}/Statistics/Mapping_and_Variant_Statistics.tab" ]; then
        mtbseq_exists=1
    fi

    """
}