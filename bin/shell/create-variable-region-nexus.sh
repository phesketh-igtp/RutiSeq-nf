#!/bin/bash
set -euo pipefail

## create-variable-region-nexus.sh
##  @author: Poppy J Hesketh Best
##  @version: v1.1.0
##  @date: 2025-04-11
##  @description:
##      This script uses the MTBSeq outputs, and a cluster file, 
##          cluster identifiers to extract the SNPs unique to the 
##          selection of genomes within a cluster and creates those
##          alignments, before extracting those genomic positions from
##          the reference genomes, and the MTBC_ancestor (Comas 2013).
##      A nexus file will be generated
##  @changelog
##         v1.0.0-2024-12-01: Initial version
##         v1.1.0-2025-04-11: Changed - Included input flags

# Help message
usage() {
    echo "Usage: $0 -c <clusterID> -p <pairwise_clusters.tsv> -f <snp_fasta.fa> -t <snp_table.tab> -m <mtbc_ancestor_path>"
    echo
    echo "Required arguments:"
    echo "  -c    Cluster ID"
    echo "  -p    Pairwise clusters file"
    echo "  -f    SNP FASTA file"
    echo "  -t    SNP annotation table"
    echo "  -m    Path to gzipped MTBC ancestor .pos.gz file"
    echo
    echo "Example:"
    echo "  $0 -c cluster123 -p pairwise.tsv -f snps.fasta -t snp_table.tab -m MTB_ancestor.pos.gz"
    exit 1
}

# Parse flags
while getopts ":c:p:f:t:m:h" opt; do
    case ${opt} in
        c) clusterID="$OPTARG" ;;
        p) pairwise_clusters="$OPTARG" ;;
        f) snp_fasta="$OPTARG" ;;
        t) snp_tab="$OPTARG" ;;
        m) mtbc_ancestor_path="$OPTARG" ;;
        h) usage ;;
        \?) echo "❌ Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "❌ Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Check required arguments
if [ -z "${clusterID:-}" ] || [ -z "${pairwise_clusters:-}" ] || [ -z "${snp_fasta:-}" ] || [ -z "${snp_tab:-}" ] || [ -z "${mtbc_ancestor_path:-}" ]; then
    echo "❌ Missing required arguments." >&2
    usage
fi

# Create output directories
mkdir -p nexus/ fasta/ positions/

# Create list of genomes
grep "${clusterID}" "${pairwise_clusters}" | cut -f1 > "${clusterID}.genomes.list"

# Extract FASTA for the cluster
> "${clusterID}.fasta"
while IFS=";" read -r genome; do
    seqkit grep -w 0 -n -p "${genome}" "${snp_fasta}" >> "${clusterID}.fasta"
done < "${clusterID}.genomes.list"

# Run snp-sites
snp-sites "${clusterID}.fasta" > "${clusterID}.snpsites.fasta"
snp-sites "${clusterID}.fasta" -v | cut -f2 | sed '1,4d' > "positions/${clusterID}_positions.tab"

# H37Rv reference
awk 'NR==FNR {pos[$1+2]; next} FNR in pos {print $3}' "positions/${clusterID}_positions.tab" "${snp_tab}" > "${clusterID}_tmp_refseq"
paste -s -d "" "${clusterID}_tmp_refseq" | sed '1i >H37Rv' > "${clusterID}_H37Rv.fasta"

# Genomic positions
awk 'NR==FNR {pos[$1+2]; next} FNR in pos {print $1}' "positions/${clusterID}_positions.tab" "${snp_tab}" > "positions/${clusterID}_genomic_positions.tab"

# MTBC ancestor
cp "${mtbc_ancestor_path}" tmp.MTB_anc.pos.gz
gunzip tmp.MTB_anc.pos.gz
awk 'NR==FNR {pos[$1]; next} FNR in pos {print $3}' "positions/${clusterID}_genomic_positions.tab" tmp.MTB_anc.pos > "${clusterID}_tmp_MTB_anc"
paste -s -d "" "${clusterID}_tmp_MTB_anc" | sed '1i >MTB_anc' > "${clusterID}_MTB_anc.fasta"
rm -rf tmp.*

# Combine FASTA files
cat "${clusterID}.snpsites.fasta" "${clusterID}_H37Rv.fasta" "${clusterID}_MTB_anc.fasta" > "fasta/${clusterID}_refseq.fasta"

# Convert to NEXUS
seqret -osformat2 nexus -sequence "fasta/${clusterID}_refseq.fasta" -outseq "nexus/${clusterID}_refseq.nex"

echo "✅ Done: NEXUS file saved to nexus/${clusterID}_refseq.nex"
