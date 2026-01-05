process SNIPPY_DATED_PHYLOGENY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        Generate timetrees for each lineage using the alignments and trees generated in the previous step.
*/

    conda params.phylogeny_env

    publishDir "${params.outDir}/db/comparison/snippy/", mode: 'copy'

    input:
        tuple path(variant_aln), 
            path(invariant_sites),
            path(contree)
        path(metadata)

    output:
        path("timetree/*"), optional: true, 
                            emit: timetrees_out

    script:

    """
    # get genome IDs from fasta
        grep '>' ${variant_aln} \\
        | sed 's@>@@g' \\
        > genomes.list

    # Isolate the sames and sampleIDs
        echo "name\tdate" > dates.tsv
        grep -f genomes.list  ${metadata} \\
                | awk -F ',' '{print \$2,\$3}' \\
                | sed 's@,@\t@g' \\
                >> dates.tsv
        sed -i 's@ @\t@g' dates.tsv
    
    # Run the timetree
        treetime \\
            --aln ${variant_aln} \\
            --tree ${contree} \\
            --dates dates.tsv \\
            --outdir timetree/ \\
            --clock-std-dev 0.2 \\
            --reroot 'least-squares'
    """

}