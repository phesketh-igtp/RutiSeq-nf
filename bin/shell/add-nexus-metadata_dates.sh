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
        m) metadata="$OPTARG" ;;
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

dates_nexus="${prefix}_dates.nexus"
locs_nexus="${prefix}_locs.nexus"

############################################################################
# Process Dates
############################################################################

cut -f2,3 -d $'\t' "$metadata" > dates.csv
n_genomes=$(wc -l < dates.csv)

Rscript -e '
d <- read.csv("dates.csv", header=FALSE)
colnames(d) <- c("sampleID", "date")
m <- table(d$sampleID, d$date)
write.table(ifelse(m > 0, 1, 0), file="out_dates_matrix.csv", sep="\t", col.names=NA, quote=FALSE)
'

sed 's/^\([^\t]*\)\t/\1 /' out_dates_matrix.csv > final_dates_matrix.out

matrix=$(cat final_dates_matrix.out)
dates_list=$(head -n 1 out_dates_matrix.csv | cut -f2- | tr '\t' ' ')

cat "$orig_nexus" > "$dates_nexus"
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
n_genomes=$(wc -l < locs.csv)

Rscript -e '
d <- read.csv("locs.csv", header=FALSE)
colnames(d) <- c("sampleID", "location")
m <- table(d$sampleID, d$location)
write.table(ifelse(m > 0, 1, 0), file="out_locs_matrix.csv", sep="\t", col.names=NA, quote=FALSE)
'

sed 's/^\([^\t]*\)\t/\1 /' out_locs_matrix.csv > final_locs_matrix.out

matrix=$(cat final_locs_matrix.out)
locs_list=$(head -n 1 out_locs_matrix.csv | cut -f2- | tr '\t' ' ')

cat "$orig_nexus" > "$locs_nexus"
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
