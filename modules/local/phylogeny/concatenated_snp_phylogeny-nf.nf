process CONCATENATED_VARIABLE_REGION_PHYLOGENY {

    tag "${runID}: ${lineage}"

    conda params.phylogeny_env

    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/", mode: 'copy'

    input:
        val(runID)
        tuple val(lineage), 
                path(fasta), 
                path(tab)

    output:
        path("Phylogeny/*")
        
        tuple val(lineage), path("Phylogeny/${lineage}_ML.contree"),
                            //path("Phylogeny/${lineage}_timetree/timetree.nexus"),
                            path("Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.aln.fasta"), emit: phylogeny_plotting_ch

    script:

    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
        mkdir Phylogeny/

        # Need to make the respective SNP alignment for the H37Rv and the Ancestral sequence for the phylogeny
            
        # 1. grab a single sequence in the fasta file (first) to get the positions
            seqkit sample -n 1 ${fasta} > ${lineage}.tmp.fasta
            
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
                cut -f1 ${lineage}.tmp.fasta_positions.tab > ${lineage}.tmp.fasta_positions
                rm ${lineage}.tmp.fasta_positions.tab
        
        # 3. obtain the reference positions (H37Rv) for the cluster positions
            for i in `cat ${lineage}.tmp.fasta_positions`; do 
                sed -n \$((i+2))'p' ${tab} | cut -f3
            done > ${lineage}.tmp_refseq
        
        # 4. convert column into fasta
            paste -s -d "" ${lineage}.tmp_refseq | sed '1i >H37Rv' > Phylogeny/${lineage}.ref-H37Rv.fasta
            
        # 5. get the genomic positions of the SNPs
            while read -r position; do
                sed -n \$((position+2))'p' ${tab} | cut -f 1; 
            done < ${lineage}.tmp.fasta_positions > Phylogeny/${lineage}_genomic_positions
            cp ${params.mtbc_ancestor_path} ${lineage}.tmp.MTB_anc.pos.gz; gunzip ${lineage}.tmp.MTB_anc.pos.gz
            
        # 6. Get the same SNPs for the 'ancestor' genomes
            for i in `cat Phylogeny/${lineage}_genomic_positions`; do 
                sed -n \${i}'p' ${lineage}.tmp.MTB_anc.pos | cut -f3 # doesnt need to +2 as the tsv file has no header
            done > ${lineage}.tmp.MTB_anc
        
        # 7. convert the column in fasta
            paste -s -d "" ${lineage}.tmp.MTB_anc | sed '1i >MTB_anc' > Phylogeny/${lineage}.ref-MTB_anc.fasta

        # 8. Merge all the sequences into a single fasta file
            cat ${fasta} \\
                Phylogeny/${lineage}.ref-H37Rv.fasta \\
                Phylogeny/${lineage}.ref-MTB_anc.fasta \\
                > Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.fasta
            
            # remove all temporary files
                rm ${lineage}.tmp.*

        # 9. Perform alignment of sequences 
            mafft --auto --thread ${params.cpu} \\
                    Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.fasta \\
                    > Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.aln.fasta

        # 10. Perform phylogeny
            iqtree -s Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.aln.fasta \\
                    -m GTR+G4 -T AUTO \\
                    -ntmax ${params.cpu} \\
                    -B ${params.iqtree_bootstraps} \\
                    --prefix ${lineage}_ML

        # Create molecular timetree
            #treetime --aln Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.aln.fasta \\
            #        --tree Phylogeny/${lineage}_ML.contree \\
            #        --dates ${params.metadata} \\
            #        --outdir Phylogeny/${lineage}_timetree

    """
}