#!/bin/bash

set -euo pipefail

# Help function
usage() {
    echo "Usage: $0 -c <clusterID> -p <pairwise_clusters> -f <snp_fasta> -t <snp_tab> -m <mtbc_ancestor_path> -n <ancestor> -l <lineage>"
    echo
    echo "Options:"
    echo "  -c    Cluster ID"
    echo "  -p    Pairwise cluster TSV (col1=sampleID, colx=clusterID)"
    echo "  -f    SNP FASTA file"
    echo "  -t    SNP table"
    echo "  -m    MTBC ancestor file path"
    echo "  -n    Ancestor file (CSV format)"
    echo "  -l    Lineage identifier"
    echo "  -h    Show this help message"
    exit 1
}

# Parse options
while getopts ":c:p:f:t:m:n:l:h" opt; do
  case ${opt} in
    c ) clusterID=$OPTARG ;;
    p ) pairwise_clusters=$OPTARG ;;
    f ) snp_fasta=$OPTARG ;;
    t ) snp_tab=$OPTARG ;;
    m ) mtbc_ancestor_path=$OPTARG ;;
    n ) ancestor=$OPTARG ;;
    l ) lineage=$OPTARG ;;
    h ) usage ;;
    \? ) echo "Invalid option: -$OPTARG" >&2; usage ;;
    : ) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Check for required arguments
if [[ -z "${clusterID:-}" || -z "${pairwise_clusters:-}" || -z "${snp_fasta:-}" || -z "${snp_tab:-}" || -z "${mtbc_ancestor_path:-}" || -z "${ancestor:-}" || -z "${lineage:-}" ]]; then
    echo "Missing required arguments." >&2
    usage
fi

# Create output and temp directories
mkdir -p nexus/ fasta/ positions/

# Create list of genomes in the cluster
grep "${clusterID}" "${pairwise_clusters}" | cut -f1 > "${clusterID}.genomes.list"

# Build cluster FASTA
> "${clusterID}.fasta"
while IFS=";" read -r genome; do
    seqkit grep -w 0 -n -p "${genome}" "${snp_fasta}" >> "${clusterID}.fasta"
done < "${clusterID}.genomes.list"

# Run snp-sites
snp-sites "${clusterID}.fasta" > "${clusterID}.snpsites.fasta"
snp-sites "${clusterID}.fasta" -v | cut -f2 | sed '1,4d' > "positions/${clusterID}_positions.tab"

# Generate NODE ancestor sequence
sed -i 's/ /,/g' "${ancestor}"
awk -F ',' 'NR==FNR {pos[$1+1]; next} FNR in pos {print $4}' "positions/${clusterID}_positions.tab" "${ancestor}" > "${clusterID}_node_anc"
paste -s -d "" "${clusterID}_node_anc" | sed "1i >MRCA" > "${clusterID}_MRCA.fasta"

# H37Rv variant positions
awk 'NR==FNR {pos[$1+2]; next} FNR in pos {print $3}' "positions/${clusterID}_positions.tab" "${snp_tab}" > "${clusterID}_tmp_refseq"
paste -s -d "" "${clusterID}_tmp_refseq" | sed '1i >H37Rv' > "${clusterID}_H37Rv.fasta"

# Genomic positions
awk 'NR==FNR {pos[$1+2]; next} FNR in pos {print $1}' "positions/${clusterID}_positions.tab" "${snp_tab}" > "positions/${clusterID}_genomic_positions.tab"

# Valencian ancestor (MTB_anc) variant positions
cp "${mtbc_ancestor_path}" "${lineage}.tmp.MTB_anc.pos.gz"
gunzip "${lineage}.tmp.MTB_anc.pos.gz"

awk 'NR==FNR {lines[$1]; next} FNR in lines {print $3}' "positions/${clusterID}_genomic_positions.tab" "${lineage}.tmp.MTB_anc.pos" > "${clusterID}_tmp_MTB_anc"
paste -s -d "" "${clusterID}_tmp_MTB_anc" | sed '1i >MTB_anc' > "${clusterID}_MTB_anc.fasta"
rm -f "${lineage}.tmp.MTB_anc.pos"

# Create final combined FASTA
cat "${clusterID}.snpsites.fasta" "${clusterID}_H37Rv.fasta" "${clusterID}_MTB_anc.fasta" "${clusterID}_MRCA.fasta" > "fasta/${clusterID}_refseq_mrca.fasta"

# Convert to NEXUS format
seqret -osformat2 nexus -sequence "fasta/${clusterID}_refseq_mrca.fasta" -outseq "nexus/${clusterID}_refseq_mrca.nex"

echo "Process complete for cluster ${clusterID}."