#!/bin/bash

##  @author: Poppy J Hesketh Best
##  @date: 2025-04-24
##  @version: 1.0.0
##  @description:
##      Providing a list of genomes, a global-SNP alignment and tabular
##          file of the 
##  @changelog:
##      v1.0.0-2025-04-15 - Initial version

############################################################################
# Parse arguments
############################################################################

while getopts "i:t:g:p:" opt; do
    case $opt in
        i) global_snp_fas="$OPTARG" ;;
        t) global_snp_tab="$OPTARG" ;;
        g) genome_list="$OPTARG" ;;
        *) echo "Usage: $0 -i global_snps.fasta -t global_snps.tab -p genome_list" >&2; exit 1 ;;
    esac
done

# Input flag check
if [ -z "$global_snp_fas" ] || [ -z "$global_snp_tab" ] || [ -z "$genome_list" ]; then
    echo "Missing arguments!"
    echo "Usage: $0 -i global_snps.fasta -t global_snps.tab -p genome_list"
    exit 1
fi

############################################################################
# Iterate over the MJN
############################################################################

# Create output directories
    mkdir -p nexus/ fasta/ positions/ raw/ tabs/

while IFS= read -r genome; do

    grep -v ${genome} ${genome_list} > ".tmp.genome.list"

    # Remove a genomes from the "global_snp_fas"
        seqkit grep -w 0 -f .tmp.genome.list ${global_snp_fas} > wo_${genome}.fasta

    # Run snp-sites
        #coresnpfilter -e -c 0.95 C1.fasta --table C1.coresnpfilt.95.pos > C1.coresnpfilt.95.fasta
        #awk '{ if ( $9 == 1 ) print $1 }' C1.coresnpfilt.95.pos > C1.coresnpfilt.95.local.pos

    # Run snp-sites
        #snp-sites "fasta/wo_${genome}.fasta" > "fasta/wo_${genome}.snpsites.fasta"
        #snp-sites "fasta/wo_${genome}.fasta" -v | cut -f2 | sed '1,4d' > "positions/wo_${genome}_positions.tab"

        coresnpfilter -e -c 1 wo_${genome}.fasta \
                --table positions/wo_${genome}.coresnpfilt.pos \
                > fasta/wo_${genome}.snpsites.fasta

        awk '{ if ( $9 == 1 ) print $1 }' positions/wo_${genome}.coresnpfilt.pos \
            > positions/wo_${genome}_positions.tab

        awk '{ if ( $9 == 1 ) print }' positions/wo_${genome}.coresnpfilt.pos \
            > raw/wo_${genome}_positions.tab

    # Genomic positions
        awk 'NR==FNR {pos[$1+2]; next} FNR in pos {print $1}' \
            "positions/wo_${genome}_positions.tab" \
            "${global_snp_tab}" > "raw/wo_${genome}_genomic_positions.tab"

        echo "genomic_pos   local_pos	a	c	g	t	count	frac	var	keep" > tabs/wo_${genome}.snp.binary.tsv
        paste raw/wo_${genome}_genomic_positions.tab raw/wo_${genome}_positions.tab > tabs/tmp.wo_${genome}.snp.binary.tsv
        cat tabs/tmp.wo_${genome}.snp.binary.tsv >> tabs/wo_${genome}.snp.binary.tsv; rm tabs/tmp.wo_${genome}.snp.binary.tsv

    rm .tmp.genome.list wo_${genome}.fasta

done < ${genome_list}

############################################################################
# Iterate over the MJN
############################################################################
