process GENERATE_TIMETREES {

    tag "${lineage}"

    conda params.phylogeny_env

    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/", mode: 'copy'

    input:
        tuple val(lineage), 
                path(contree), 
                path(alignments)
        path(metadata)

    output:
        path("${lineage}_timetree/*")

        tuple val(lineage), 
                path("${lineage}_timetree/timetree.nexus"),
                path("${lineage}_timetree/ancestral_sequences.fasta"), emit: timetrees_ch

    script:
    """
    # get genome IDs from fasta
        grep '>' ${alignments} | sed 's@>@@g' > genomes.list

    # Isolate the sames and sampleIDs
        grep -f genomes.list  ${metadata} \\
                | awk -F ',' '{print $2,$3}' \\
                | sed 's@,@\t@g' > dates.csv
                
    # Run the timetree
        treetime --aln ${alignments} \\
                --tree ${contree} \\
                --dates dates.csv \\
                --outdir ${lineage}_timetree
    """
}