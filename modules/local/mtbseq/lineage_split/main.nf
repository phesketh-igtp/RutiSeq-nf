process MTBSEQ_LINEAGE_SPLITTING {

    tag "${runID}_${lineage}"

    publishDir "${params.outdir}/bbdd/mtbseq/lineages/${runID}", mode: 'copy'

    input:
    val runID
    path tbprofile_tbdb_compile
    each lineage

    output:
    tuple val(runID), val(lineage), path("${lineage}.samples.txt"), emit: lineages_ch
    path "${runID}.unique-lineages.list",                           emit: unique_lineages

    script:
    """
    awk '{if ($3 == "${lineage}") print}'  ${tbprofile_tbdb_compile} | cut -f1 > ${lineage}.samples.txt

    # Create a list of unique lineages found in this run
    if [ ! -f "${runID}.unique-lineages.list" ]; then
        echo "${lineage}" > "${runID}.unique-lineages.list"
    else
        echo "${lineage}" >> "${runID}.unique-lineages.list"
    fi

    # Sort and remove duplicates from the unique lineages list
    sort -u -o "${runID}.unique-lineages.list" "${runID}.unique-lineages.list"
    
    """

}