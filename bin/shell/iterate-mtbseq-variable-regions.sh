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
    mkdir -p Amend/ Joint/ Groups/ Called/ Position_Tables Results/

    while IFS= read -r genome; do

        grep -v ${genome} ${genome_list} > ".tmp.genome.list"
        mkdir -p Results/wo_${genome}/

        for sampleID in $(cat .tmp.genome.list); do
            ln -s /imppc/labs/emlab/share/TBSEQ/RutiSeq/bbdd/mtbseq/samples/${sampleID}/Called/* Called/
            ln -s /imppc/labs/emlab/share/TBSEQ/RutiSeq/bbdd/mtbseq/samples/${sampleID}/Position_Tables/* Position_Tables/
        done

        MTBseq --step TBJoin --thread 8 --minbqual 20 --mincovf 4 \
                --mincovr 4 --minphred20 4 --minfreq 75 --unambig 95 --window 10
        MTBseq --step TBAmend --thread 8 --minbqual 20 --mincovf 4 \
                --mincovr 4 --minphred20 4 --minfreq 75 --unambig 95 --window 10

        cat Amend/*_w10.fasta > Results/wo_${genome}/wo_${genome}.fasta
        cat Amend/*_w10.tab >   Results/wo_${genome}/wo_${genome}.tab

        # Run snp-sites
        snp-sites "Results/wo_${genome}/wo_${genome}.fasta" > "fasta/wo_${genome}.snpsites.fasta"
        snp-sites "Results/wo_${genome}/wo_${genome}.fasta" -v \
            | cut -f2 | sed '1,4d' > "Results/wo_${genome}/wo_${genome}_positions.tab"

        # Genomic positions
        awk 'NR==FNR {pos[$1+2]; next} FNR in pos {print $1}' \
            Results/wo_${genome}/wo_${genome}_positions.tab \
            Results/wo_${genome}/wo_${genome}.tab \
            > Results/wo_${genome}/wo_${genome}_genomic_positions.tab

        paste Results/wo_${genome}/wo_${genome}_positions.tab \
            Results/wo_${genome}/wo_${genome}_genomic_positions.tab \
            > Results/wo_${genome}/wo_${genome}_positions.tab

        rm Results/wo_${genome}/.tmp.genome.list

done < ${genome_list}

############################################################################
# Iterate over the MJN
############################################################################
