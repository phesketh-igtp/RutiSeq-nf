process MTBSEQ_LINEAGE_SPLITTING {

    tag "${runID}_${lineage}"

    publishDir "${params.outdir}/bbdd/mtbseq/lineages/${runID}", mode: 'copy'

    input:
    val runID
    path tbprofile_tbdb_compile
    each lineage

    output:
    tuple val(runID), val(lineage), path("${lineage}.samples.txt"),         emit: lineages_ch
    path "${runID}.unique-lineages.list",                                   emit: unique_lineages
    path "samples-mtbseq.csv"
    tuple val(lineage), val(sample_name), path(called_dir), path(position_tables_dir),    emit: mtbseq_paths

    script:
    """
    # Create the lineage-specific sample list
    awk '{if (\$3 == "${lineage}") print \$1}'  ${tbprofile_tbdb_compile} > ${lineage}.samples.txt

    # Create samples-mtbseq.csv
    paste -d ',' ${lineage}.samples.txt samples.called samples.pos > samples-mtbseq.csv

    # Create the tuple for mtbseq_paths_ch
    while IFS=',' read -r sample_name called_dir position_tables_dir; do
        echo "${lineage},\$sample_name,\$called_dir,\$position_tables_dir"
    done < samples-mtbseq.csv > mtbseq_paths.txt

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