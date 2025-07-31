#!/bin/bash

##  @author: Poppy J Hesketh Best
##  @date: 2025-04-15
##  @version: 1.0.0
##  @description:
##  @changelog:
##      v1.0.0-2025-04-15 - Initial version

############################################################################
# Parse arguments
############################################################################

while getopts "i:m:p:" opt; do
    case $opt in
        i) orig_nexus="$OPTARG" ;;
        p) prefix="$OPTARG" ;;
        *) echo "Usage: $0 -i input.nexus -m metadata.tsv -p output_prefix" >&2; exit 1 ;;
    esac
done

# Sanity check
if [ -z "$orig_nexus" ] || [ -z "$metadata" ] || [ -z "$prefix" ]; then
    echo "Missing arguments!"
    echo "Usage: $0 -i input.nexus -m metadata.tsv -p output_prefix"
    exit 1
fi

dates_nexus="${prefix}_dates.nex"
locs_nexus="${prefix}_locs.nex"

############################################################################
# Process Dates
############################################################################

N_SAMPLES=$(( 2 + $(wc -l < dates.csv) ))

Rscript -e '
    d <- read.csv("dates.csv", header=FALSE)
    colnames(d) <- c("sampleID", "date")

    d$date <- format(as.Date(d$date), "%Y")

    m <- table(d$sampleID, d$date)

    write.table(ifelse(m > 0, 1, 0), 
        file="out_dates_matrix.csv", 
        sep=",", col.names=NA, quote=FALSE)
'

dates_list=$( sed 's@^,@@g' out_dates_matrix.csv | head -1 | sed 's@,@ @g' )
matrix=$( sed 's/,/\t/' out_dates_matrix.csv | sed '1d' | cat )

cat "$orig_nexus" > "$dates_nexus"
sed -i '/begin assumptions;/,/end;/d' "$dates_nexus"

echo -e "
BEGIN TRAITS;
Dimensions NTRAITS=${n_genomes};
Format labels=yes missing=? separator=Comma;
TraitLabels ${dates_list};

Matrix
${matrix}
;

END;
" >> "$dates_nexus"

######################################################################

TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
N_SAMPLES= # numeric value, e.g. 5
N_NUCLEOTIDES= # numeric value, e.g. 5
ALN_DATA_MATRIX= # must be tab separated matrix e.g. : SAMPLEID ACGTACGTACGT
N_YEARS= # numeric value, e.g. 5
YEAR_LIST= # must be whitespace seperated string e.g. : 2021 2022 2023 2024
DATE_DATA_MATRIX=

######################################################################

echo -e "#NEXUS
[TITLE: Outbreak cluster nexus with metadata, ${TIMESTAMP}]

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
    Dimensions NTRAITS=${N_YEARS};
    Format labels=yes missing=? separator=Comma;
    TraitLabels ${YEAR_LIST};
MATRIX

${DATE_DATA_MATRIX}
;
END
" > ${prefix_out}.nex
