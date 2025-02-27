process CONCATENATED_VARIABLE_REGION_PHYLOGENY {

    tag "${runID}: ${lineage}"

    conda params.phylogeny_env

    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/", mode: 'copy'

    input:
        val(runID)
        tuple val(lineage), 
                path(fasta), 
                path(tab)

    output:
        path("Phylogeny/*")
        
        tuple val(lineage), path("Phylogeny/${lineage}_ML.contree"),
                            path("Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.aln.fasta"), emit: phylogeny_plotting_ch

    script:

    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    # Create the fasta files for the phylogeny
        bash ${params.script_dir}/shell/concatenate-variable-pylogeny-ancestors.sh \\
                ${fasta} ${lineage} \\
                ${tab} ${params.mtbc_ancestor_path}

        # Perform alignment of sequences 
            mafft --auto --thread ${params.cpus} \\
                    Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.fasta \\
                    > Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.aln.fasta

        # Perform phylogeny
            iqtree -s Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.aln.fasta \\
                    -m ${params.iqtree_model} \\
                    -T AUTO \\
                    -ntmax ${params.cpus} \\
                    -B ${params.iqtree_bootstraps} \\
                    --prefix ${lineage}_ML

            mv ${lineage}_ML.* Phylogeny/
    """

}