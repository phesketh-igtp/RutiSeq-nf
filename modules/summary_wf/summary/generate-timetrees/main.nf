process DATED_PHYLOGENY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        Generate timetrees for each lineage using the alignments and trees generated in the previous step.
*/

    tag "${lineage}"

    conda params.phylogeny_env

    publishDir "${params.outDir}/db/mtbseq/pairwise/${lineage}/", mode: 'copy'

    input:
        tuple val(lineage), 
                path(fasta), 
                path(tab)
        path(metadata)

    output:
        path("${lineage}_timetree/")

        tuple val(lineage), 
                path("${lineage}_timetree/timetree.nexus"),
                path("${lineage}_timetree/ancestral_sequences.fasta"), 
                            optional: true, 
                            emit: timetrees_out

    script:
    """
    # get genome IDs from fasta
        grep '>' ${fasta} | sed 's@>@@g' > genomes.list

    # Isolate the sames and sampleIDs
        echo "name\tdate" > dates.tsv
        
        grep -f genomes.list  ${metadata} \\
                | awk -F ',' '{print \$2,\$3}' \\
                | sed 's@,@\t@g' \\
                >> dates.tsv

        sed -i 's@ @\t@g' dates.tsv

    # Perform main phylogeny
        iqtree -s ${fasta} \\
        -m ${params.iqtree_model} \\
        -T AUTO \\
        -ntmax ${params.cpus} \\
        -B ${params.iqtree_bootstraps} \\
        --prefix ${lineage}_reference-free

    mkdir ${lineage}_timetree/

    # Run the timetree
        treetime \\
            --aln ${fasta} \\
            --tree ${lineage}_reference-free.contree \\
            --dates dates.tsv \\
            --outdir ${lineage}_timetree/ \\
            --clock-std-dev 0.2 \\
            --reroot 'least-squares' || \\
    echo "TreeTime failed for ${lineage}, skipping time tree generation" >&2
    """

}