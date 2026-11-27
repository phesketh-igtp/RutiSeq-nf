#!/usr/bin/env bash
#$ -N treetime
#$ -cwd
#$ -j y
#$ -o treetime.$JOB_ID.log
#$ -e treetime.$JOB_ID.err
#$ -q d12imppc
#$ -pe smp 8
#$ -l h_rt=16:00:00
# ----------------------------
# History
# ----------------------------
#' Author: Poppy J Hesketh Best
#' Date: 2026-05-08
#' Changelog:
#' - v1.0: Initial version (2026-05-08)
version="1.0 (2026-05-08)"
last_modified_date="2026-05-08"

source /imppc/labs/emlab/phesketh/miniconda3/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate gubbins

###############################################################
# Link all the VCF files
echo "Log: Collecting all VCF files..."
###############################################################
# Get updated list of samples in DB
outDir="/imppc/labs/emlab/share/TBSEQ/RutiSeq"
ls ${outDir}/db/samples/ | sort | uniq > samples.list

# Assign the base path
BASE="${outDir}/db/samples"

mkdir -p src/

while read -r f; do
    src_dir="\$BASE/\$f/snippy-ng/"

    for file in "\$src_dir"/*.QC.mask.filt.vcf.gz \\
                "\$src_dir"/*.QC.mask.filt.snps.consensus.fa.gz; do

        # Only create symlink if file exists
        if [[ -f "\$file" ]]; then
            ln -sf "\$file" src/
            echo -e "\$f" >> samples.list.new
        fi

    done
done < <(sort -u samples.list)

###############################################################
# Concatenate the consensus sequences
echo "Log: Concatenating consensus sequence..."
###############################################################
cat src/*.QC.mask.filt.snps.consensus.fa.gz \\
    > RutiSeq.QC.mask.filt.snps.consensus.fa.gz

gunzip --keep RutiSeq.QC.mask.filt.snps.consensus.fa.gz

seqkit rmdup -n RutiSeq.QC.mask.filt.snps.consensus.fa \\
    >> RutiSeq.QC.mask.filt.snps.consensus.fa.tmp
mv RutiSeq.QC.mask.filt.snps.consensus.fa.tmp \\
    RutiSeq.QC.mask.filt.snps.consensus.fa

###############################################################
# Run gubbins
echo "Log: Performing gubbins analysis..."
###############################################################
run_gubbins.py \\
    --extensive-search \\
    --prefix RutiSeq.gubbins \\
    --threads ${task.cpus} \\
    --verbose \\
    RutiSeq.QC.mask.filt.snps.consensus.fa

gff2bed < RutiSeq.gubbins.recombination_predictions.gff \\
    > RutiSeq.gubbins.recombination_predictions.bed

sed -i 's@SEQUENCE@NC_000962.3@g' \\
    RutiSeq.gubbins.recombination_predictions.bed
    
# Clean up
mkdir gubbins/; mv RutiSeq.gubbins.* gubbins/

###############################################################
# Filter VCFs with gubbins masking sites and regen consensus
echo "Log: Generating gubbins filtering VCF file..."
###############################################################
for vcf in src/*.QC.mask.filt.vcf.gz; do

    src_dir="${outDir}/${f}/snippy-ng/"
    f=$(basename $f .QC.mask.filt.vcf.gz)

    bcftools filter \\
        src/\${f}.QC.mask.filt.vcf.gz \\
        --mask-file gubbins/RutiSeq.gubbins.recombination_predictions.gff \\
        --soft-filter "gubbins" \\
        -Oz -o src/\${f}.QC.mask.filt.gubbins.core.vcf.gz

    # Apply only SNPs to consensus genome
    bcftools view \\
        -f PASS \\
        -v snps \\
        src/\${f}.QC.mask.filt.gubbins.core.vcf.gz \\
        -Oz \\
        -o src/\${f}.QC.mask.filt.gubbins.snps.vcf.gz

    tabix -f src/\${f}.QC.mask.filt.gubbins.snps.vcf.gz

    # Concensus
    bcftools consensus \\
        -f ${snippy_reference_fa} \\
        src/\${f}.QC.mask.filt.gubbins.snps.vcf.gz \\
        > src/\${f}.QC.mask.filt.gubbins.snps.core.consensus.fa

    sed -i "s@>NC_000962.3@>\${f}@g" \\
        src/\${f}.QC.mask.filt.gubbins.snps.core.consensus.fa

    gzip --best --keep src/\${f}.QC.mask.filt.gubbins.snps.core.consensus.fa

    # Copy over new vcfs and consensus files to sample directory
    cp src/\${f}.QC.mask.filt.gubbins.core.vcf.gz \${src_dir}/
    cp src/\${f}.QC.mask.filt.gubbins.snps.core.consensus.fa.gz \${src_dir}/

done

###############################################################
# Concatenate the consensus genomes
echo "Log: Concatenate gubbins consensus sequence..."
###############################################################
cat src/*.QC.mask.filt.gubbins.snps.consensus.fa \\
    > RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.fa

gzip --best --keep RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.fa

###############################################################
# Calculate the pairwise distances between the consensus genomes
echo "Log: Calculating the pairwise distance of all consensus sequences..."
###############################################################

pdsm --threads ${task.cpus} --sort --progress \\
    RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.fa.gz \\
    > RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.pdsm.csv

pdsm --threads ${task.cpus} --sort --progress \\
    RutiSeq.QC.mask.filt.snps.consensus.fa \\
    > RutiSeq.QC.mask.filt.snps.consensus.pdsm.csv
