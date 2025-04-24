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

n_genomes=$(( 2 + $(wc -l < dates.csv) ))

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

############################################################################
# Process Locations
############################################################################

cut -f2,4 -d $'\t' "$metadata" > locs.csv

Rscript -e '
    d <- read.csv("locs.csv", header=FALSE)

    colnames(d) <- c("sampleID", "location")

    m <- table(d$sampleID, d$location)

    write.table(ifelse(m > 0, 1, 0), 
        file="out_locs_matrix.csv", 
        sep=",", col.names=NA, quote=FALSE)
'

locs_list=$( sed 's@^,@@g' out_locs_matrix.csv | head -1 | sed 's@,@ @g' )
matrix=$( sed 's/,/\t/' out_locs_matrix.csv | sed '1d' | cat )

cat "$orig_nexus" > "$dates_nexus"
sed -i '/begin assumptions;/,/end;/d' "$dates_nexus"

echo -e "
BEGIN TRAITS;
Dimensions NTRAITS=${n_genomes};
Format labels=yes missing=? separator=Comma;
TraitLabels ${locs_list};

Matrix
${matrix}
;

END;
" >> "$locs_nexus"
