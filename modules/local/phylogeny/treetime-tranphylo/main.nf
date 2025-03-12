process TREETIME_TRANSPHYLO_ANALYSIS {

    tag "${clusterID}: ${lineage}"

    conda params.phylogeny_env

    publishDir "${params.outdir}/results/${runID}/TransPhylo/", mode: 'copy'

    input:
        val(runID)
        tuple val(lineage),
                val(clusterID),
                path(fasta), 
                path(tab)
        path(processed_clusters)
        path(metadata)

    output:
        path("${clusterID}/*")

    script:

    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    # select the clsuter genomes and filter to retain their variant sites FASTA
    grep ${clusterID} ${processed_clusters} | cut -f1 > tmp.cluster.genomes

    while IFS=";" read -r genome; do
        seqkit grep -w 0 -n -p \${genome} ${fasta} >> ${clusterID}.fasta
    done < tmp.cluster.genomes

    # Perform alignment of sequences 
        mafft --auto --thread ${params.cpus} \\
                ${clusterID}.fasta \\
                > ${clusterID}.aln.fasta

    # Perform phylogeny
        iqtree -s ${clusterID}.aln.fasta \\
                -m ${params.iqtree_model} \\
                -T AUTO \\
                -ntmax ${params.cpus} \\
                -B ${params.iqtree_bootstraps} \\
                --prefix ${clusterID}_ML

    # Isolate the sames and sampleIDs
        echo "name\tdate" > dates.tsv
        
        grep -f tmp.cluster.genomes ${metadata} \\
                | awk -F ',' '{print \$2,\$3}' \\
                | sed 's@,@\t@g' >> dates.tsv

        sed -i 's@ @\t@g' dates.tsv

    # Run the timetree
        treetime --aln ${clusterID}.aln.fasta \\
                --tree ${clusterID}_ML.contree \\
                --dates dates.tsv \\
                --outdir ${clusterID}_timetree

    mkdir -p ${clusterID}/; rm *tmp*
    mv ${clusterID}* ${clusterID}/

    """

}