#!/bin/bash

##  @author: Poppy J Hesketh Best
##  @date: 2025-07-31
##  @version: 1.0.0
##  @description: This script processes a Nexus file to add meta about sample locs. 
##                 It will produce a new Nexus file with the locs included in a 
##                 traits section.
##  @changelog:
##      v1.0.0-2025-07-31 - Initial version (draft by Poppy J Hesketh Best)

############################################################################
# Parse arguments
############################################################################

while getopts "i:m:p:" opt; do
    case $opt in
        i) aln="$OPTARG" ;;
        p) prefix="$OPTARG" ;;
        m) meta="$OPTARG" ;;
        *) echo "Usage: $0 -i alignment.fasta -m meta.tsv -p prefix
        meta: two columns, (1) sampleID and (2) date (YYYY-MM-DD), should not have a header." >&2; exit 1 ;;
    esac
done

# Sanity check
if [ -z "$aln" ] || [ -z "$meta" ] || [ -z "$prefix" ]; then
    echo "Missing arguments!"
    echo "Usage: $0 -i alignment.fasta -m meta.tsv -p prefix"
    exit 1
fi


############################################################################
# Process locs
############################################################################

grep '>' ${aln} | sed 's/>//g' > tmp.sampleid.list

cut -f1 -d ',' ${meta} | sed '1d' > tmp.locs_mat.1.csv

awk -F',' 'NR > 1 {
  for (i = 2; i <= NF; i++) {
    printf "%s", $i
    if (i < NF) printf " "
  }
  print ""
}' ${meta} > tmp.locs_mat.2.csv

paste -d '\t' tmp.locs_mat.1.csv tmp.locs_mat.2.csv | \
    grep -f tmp.sampleid.list - > tmp.locs_mat.final.csv

######################################################################
## Set all the variables needed for the nexus file
######################################################################

TIMESTAMP=$( date +"%Y-%m-%d %H:%M:%S" )

SAMPLEID_LIST=$( cat tmp.sampleid.list )

N_SAMPLES=$( grep -c '>' ${aln} | cut -f1 | sed 's/ //g' )

N_NUCLEOTIDES=$( awk '!/^>/ { len = length($0); if (len > max) max = len } END { print max }' ${aln} | sed 's/ //g' )

ALN_DATA_MATRIX=$( awk '/^>/ {if (seq) print name "\t" seq; name=substr($0,2); seq=""} !/^>/ {seq=seq $0} END {print name "\t" seq}' ${aln} | sed 's/ //g' ) 

N_LOCATIONS=$( sort tmp.locs_mat.2.csv | uniq | wc -l | sed 's/ //g' )

LOCATION_MATRIX=$( cat tmp.locs_mat.final.csv )

######################################################################
# Create Nexus file with metadata
######################################################################

echo -e "#NEXUS
[TITLE: Outbreak cluster nexus with meta, ${TIMESTAMP}]

BEGIN TAXA;
    DIMENSIONS NTAX=${N_SAMPLES};
TAXLABLES
${SAMPLEID_LIST}
;

BEGIN CHARACTERS;
    DIMENSIONS NCHAR=${N_NUCLEOTIDES};
    FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=. ;
MATRIX
${ALN_DATA_MATRIX}
;
END;

BEGIN GeoTags;
    Dimensions NCLUSTS=${N_LOCATIONS};
    Format labels=yes separator=Spaces;
MATRIX
${LOCATION_MATRIX}
;
END
" > ${prefix}.annotated.nex

######################################################################
# Clean up temporary files
######################################################################

rm tmp.*