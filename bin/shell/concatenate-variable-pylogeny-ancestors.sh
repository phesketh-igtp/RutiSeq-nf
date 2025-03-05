#!/bin/bash
##set -euo pipefail

fasta=$1
lineage=$2
tab=$3
mtbc_ancestor_path=$4

mkdir -p Phylogeny/

# Need to make the respective SNP alignment for the H37Rv and the Ancestral sequence for the phylogeny

# 1. grab a single sequence in the fasta file (first) to get the positions
    seqkit seq -w 0 ${fasta} | head -2 > ${lineage}.tmp.fasta

# 2. create list of how many positions there are in the seq.)
    awk '
        BEGIN {position = 0}
            /^>/ {next}  # Skip header lines
            {
                # Process sequence lines
                for (i = 1; i <= length($0); i++) {
                    position++
                    print position "\t" substr($0, i, 1) >> "'"${lineage}.tmp.fasta_positions.tab"'"
                }
            }' "${lineage}.tmp.fasta"
                
        cut -f1 ${lineage}.tmp.fasta_positions.tab > ${lineage}.tmp.fasta_positions
        #rm ${lineage}.tmp.fasta_positions.tab

# 3. obtain the reference positions (H37Rv) for the cluster positions
    #for i in `cat ${lineage}.tmp.fasta_positions`; do 
    #    sed -n $((i+2))'p' ${tab} | cut -f3
    #done > ${lineage}.tmp_refseq

    # newer faster version
    awk 'NR==FNR {pos[$1+2]; next} FNR in pos {print $3}' ${lineage}.tmp.fasta_positions ${tab} > ${lineage}.tmp_refseq

        
# 4. convert column into fasta
    paste -s -d "" ${lineage}.tmp_refseq | sed '1i >H37Rv' > Phylogeny/${lineage}.ref-H37Rv.fasta

# 5. get the genomic positions of the SNPs
    #while read -r position; do
    #    sed -n $((position+2))'p' ${tab} | cut -f 1; 
    #done < ${lineage}.tmp.fasta_positions > Phylogeny/${lineage}_genomic_positions.tab

    # newer faster version
    awk 'NR==FNR {pos[$1+2]; next} FNR in pos {print $1}' ${lineage}.tmp.fasta_positions ${tab} > Phylogeny/${lineage}_genomic_positions.tab


    cp ${mtbc_ancestor_path} ${lineage}.tmp.MTB_anc.pos.gz; gunzip ${lineage}.tmp.MTB_anc.pos.gz

# 6. Get the same SNPs for the 'ancestor' genomes
    #for i in `cat Phylogeny/${lineage}_genomic_positions.tab`; do 
    #    sed -n ${i}'p' ${lineage}.tmp.MTB_anc.pos | cut -f3 # doesnt need to +2 as the tsv file has no header
    #done > ${lineage}.tmp.MTB_anc

    # newer faster version
    awk 'NR==FNR {pos[$1]; next} FNR in pos {print $3}' Phylogeny/${lineage}_genomic_positions.tab ${lineage}.tmp.MTB_anc.pos > ${lineage}.tmp.MTB_anc


# 7. convert the column in fasta
    paste -s -d "" ${lineage}.tmp.MTB_anc | sed '1i >MTB_anc' > Phylogeny/${lineage}.ref-MTB_anc.fasta

# 8. Merge all the sequences into a single fasta file
    cat ${fasta} Phylogeny/${lineage}.ref-H37Rv.fasta Phylogeny/${lineage}.ref-MTB_anc.fasta > Phylogeny/${lineage}.ref-H37Rv_MTBc-anc.fasta
