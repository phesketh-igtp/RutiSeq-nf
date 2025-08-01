#!/bin/bash

##  @author: Poppy J Hesketh Best
##  @date: 2025-07-31
##  @version: 1.0.0
##  @description: This script processes a Nexus file to add meta about sample dates. 
##                 It will produce a new Nexus file with the dates included in a 
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
# Process Dates
############################################################################

Rscript -e '
    meta_file <- commandArgs(TRUE)[1]
    d <- read.csv(meta_file, header=FALSE)
    colnames(d) <- c("sampleID", "date")

    d$date <- format(as.Date(d$date), "%Y")
    m <- table(d$sampleID, d$date)

    # Convert table to matrix
    m_mat <- as.matrix(m)

    # Add the "Anc" column to m_mat and initialize with 0
    m_mat <- cbind(m_mat, Anc = 0)

    # Create new rows for H37Rv and MTB_anc with all zeros
    ref <- matrix(0, nrow = 2, ncol = ncol(m_mat),
        dimnames = list(c("H37Rv", "MTB_anc"), colnames(m_mat)))

    # Set Anc = 1 for MTB_anc
    ref["MTB_anc", "Anc"] <- 1

    # Bind rows
    m_new <- rbind(m_mat, ref)

    # Convert counts to binary presence/absence (1/0)
    m_new_bin <- ifelse(m_new > 0, 1, 0)

    # Write to file
    write.table(m_new_bin,
                file = "tmp.dates_mat.csv",
                sep = ",", col.names = NA, quote = FALSE)
' "${meta}"

sed -i '' '/^$/d' tmp.dates_mat.csv

grep '>' ${aln} | sed 's/>//g' > tmp.sampleid.list

cut -f1 -d ',' tmp.dates_mat.csv | sed '1d' > tmp.dates_mat.1.csv

awk -F',' 'NR > 1 {
  for (i = 2; i <= NF; i++) {
    printf "%s", $i
    if (i < NF) printf ","
  }
  print ""
}' tmp.dates_mat.csv > tmp.dates_mat.2.csv

paste -d '\t' tmp.dates_mat.1.csv tmp.dates_mat.2.csv | \
    grep -f tmp.sampleid.list - > tmp.dates_mat.final.csv

head -1 tmp.dates_mat.csv | sed 's/,/ /g' > tmp.dates.list

######################################################################
## Set all the variables needed for the nexus file
######################################################################

TIMESTAMP=$( date +"%Y-%m-%d %H:%M:%S" )

SAMPLEID_LIST=$( cat tmp.sampleid.list )

N_SAMPLES=$( grep -c '>' ${aln} | cut -f1 | sed 's/ //g' )

N_NUCLEOTIDES=$( awk '!/^>/ { len = length($0); if (len > max) max = len } END { print max }' ${aln} | sed 's/ //g' )

ALN_DATA_MATRIX=$( awk '/^>/ {if (seq) print name "\t" seq; name=substr($0,2); seq=""} !/^>/ {seq=seq $0} END {print name "\t" seq}' ${aln} | sed 's/ //g' ) 

N_DATES=$( sed 's/ /\n/g' tmp.dates.list | sed '1d' | wc -l | sed 's/ //g' ) # numeric value, e.g. 5

DATES_LIST=$( sed 's@^,@@g' tmp.dates_mat.csv | head -1 | sed 's@,@ @g' )

DATES_DATA_MATRIX=$( cat tmp.dates_mat.final.csv )

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

BEGIN TRAITS;
    Dimensions NTRAITS=${N_DATES};
    Format labels=yes missing=? separator=Comma;
    TraitLabels ${DATES_LIST};
MATRIX

${DATES_DATA_MATRIX}
;
END
" > ${prefix}.annotated.nex

######################################################################
# Clean up temporary files
######################################################################

rm tmp.*