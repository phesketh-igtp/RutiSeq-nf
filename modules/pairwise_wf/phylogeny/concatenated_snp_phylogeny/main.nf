process SNP_PHYLOGENY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process performs the concatenated variable region phylogeny analysis
        for concatenated SNP alignments outputs from the MTBSeq pipeline. 
        It uses the MAFFT and IQ-Tree software to perform the phylogeny analysis.
*/

    tag "${lineage}"

    conda params.phylogeny_env

    publishDir "${params.outDir}/db/comparison/mtbseq/${lineage}/", mode: 'copy'

    input:
        tuple val(lineage), 
                path(fasta), 
                path(tab)

    output:
        path("Phylogeny/*")
        
        tuple val(lineage), path("Phylogeny/${lineage}_ML.contree"),
                            path("Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.aln.fasta"), 
                            emit: phylogeny_plotting_ch

    script:

    """
    mkdir -p Phylogeny/

    #--------------------------------------------------------------------------#
    # Create respective SNP alignment for the H37Rv and the Ancestral 
    ##  sequence for the phylogeny
    #--------------------------------------------------------------------------#

    # 1. grab a single sequence in the fasta file (first) to get the positions
        seqkit seq -w 0 ${fasta} \\
            | head -2 \\
            > ${lineage}.tmp.fasta

    # 2. create list of how many positions there are in the seq.)
        awk '
            BEGIN {position = 0}
                /^>/ {next}  # Skip header lines
                {
                    # Process sequence lines
                    for (i = 1; i <= length(\$0); i++) {
                        position++
                        print position "\t" substr(\$0, i, 1) >> "'"${lineage}.tmp.fasta_positions.tab"'"
                    }
                }' "${lineage}.tmp.fasta"
                    
            cut -f1 ${lineage}.tmp.fasta_positions.tab \\
                > ${lineage}.tmp.fasta_positions

    # 3. obtain the reference positions (H37Rv) for the cluster positions
        awk 'NR==FNR {pos[\$1+2]; next} FNR in pos {print \$3}' \\
            ${lineage}.tmp.fasta_positions ${tab} \\
            > ${lineage}.tmp_refseq

    # 4. convert column into fasta
        paste -s -d "" ${lineage}.tmp_refseq \\
        | sed '1i >H37Rv' \\
        > Phylogeny/${lineage}.ref-H37Rv.fasta

    # 5. get the genomic positions of the SNPs
        awk 'NR==FNR {pos[\$1+2]; next} FNR in pos {print \$1}' \\
            ${lineage}.tmp.fasta_positions ${tab} \\
            > Phylogeny/${lineage}_genomic_positions.tab

        cp ${params.mtbc_ancestor_path} \\
            ${lineage}.tmp.MTB_anc.pos.gz
        
        gunzip ${lineage}.tmp.MTB_anc.pos.gz

    # 6. Get the same SNPs for the 'ancestor' genomes
        awk 'NR==FNR {pos[\$1]; next} FNR in pos {print \$3}' \\
            Phylogeny/${lineage}_genomic_positions.tab \\
            ${lineage}.tmp.MTB_anc.pos \\
            > ${lineage}.tmp.MTB_anc

    # 7. convert the column in fasta
        paste -s -d "" ${lineage}.tmp.MTB_anc \\
            | sed '1i >MTB_anc' \\
            > Phylogeny/${lineage}.ref-MTB_anc.fasta

    # 8. Merge all the sequences into a single fasta file
        cat ${fasta} Phylogeny/${lineage}.ref-H37Rv.fasta \\
            Phylogeny/${lineage}.ref-MTB_anc.fasta \\
            > Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.fasta

    #--------------------------------------------------------------------------#
    ## Perform the alignments and phylogeny
    #--------------------------------------------------------------------------#

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

    # Move the outputs to the phylogeny directory
        mv ${lineage}_ML.* Phylogeny/
    """

}