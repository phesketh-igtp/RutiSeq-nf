process SNIPPY_LINEAGE_CORE {

    tag "${lineage}"

    conda params.gubbins_env
    
    publishDir "${params.outDir}/db/comparison/snippy/${lineage}/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple val(lineage), 
            val(sampleIDs), // Collection of VCF files from snippy runs
            val(sampleID_count)

    output:
        path("gubbins/")
        path("${lineage}.QC.mask.filt.snps.consensus.fa.gz")
        path("${lineage}.QC.mask.filt.snps.consensus.pdsm.csv")
        path("${lineage}.QC.mask.filt.gubbins.snps.consensus.fa.gz")
        path("${lineage}.QC.mask.filt.gubbins.snps.consensus.pdsm.csv")
        
        tuple val(lineage),
            path("${lineage}.QC.mask.filt.gubbins.snps.consensus.fa.gz"),
                emit: snippy_lineage_phylo_ch

    script:

    """
    ###############################################################
    # Link all the VCF files
    ###############################################################
    echo "${sampleIDs.join('\n')}" | sort | uniq > samples.list

    # Assign the base path
    BASE="${params.outDir}/db/samples"

    while read -r f; do

        src_dir="\$BASE/\$f/snippy-ng/"
        dest_dir="src/"

        # Only create dest directory if needed
        mkdir -p src/
        ln -sf "\$src_dir"/*.QC.mask.filt.vcf.gz src/
        ln -sf "\$src_dir"/*.QC.mask.filt.snps.consensus.fa.gz

    done < samples.list

    ###############################################################
    # Concatenate the consensus sequences
    ###############################################################
    cat *.QC.mask.filt.snps.consensus.fa.gz \\
        > ${lineage}.QC.mask.filt.snps.consensus.fa.gz

    gunzip --keep ${lineage}.QC.mask.filt.snps.consensus.fa.gz

    ###############################################################
    # Run gubbins
    ###############################################################
    run_gubbins.py \\
        --extensive-search \\
        --prefix ${lineage}.gubbins \\
        --threads ${task.cpus} \\
        ${lineage}.QC.mask.filt.snps.consensus.fa

    gff2bed < ${lineage}.gubbins.recombination_predictions.gff \\
        > ${lineage}.gubbins.recombination_predictions.bed

    sed -i 's@SEQUENCE@NC_000962.3@g' \\
        ${lineage}.gubbins.recombination_predictions.bed
    
    # Clean up
    mkdir gubbins/; mv ${lineage}.gubbins.* gubbins/

    ###############################################################
    # Filter VCFs with gubbins masking sites and regen consensus
    ###############################################################
    # Mask with gubbins
    while read -r f; do

        src_dir="\${BASE}/\${f}/snippy-ng/"

        bcftools filter \\
            src/\${f}.QC.mask.filt.vcf.gz \\
            --mask-file gubbins/RutiSeq.gubbins.recombination_predictions.gff \\
            --soft-filter "gubbins" \\
            -Oz -o src/\${f}.QC.mask.filt.gubbins.lineage.vcf.gz

        # Apply only SNPs to consensus genome
        bcftools view \\
            -f PASS \\
            -v snps \\
            src/\${f}.QC.mask.filt.gubbins.lineage.vcf.gz \\
            -Oz \\
            -o src/\${f}.QC.mask.filt.gubbins.snps.vcf.gz

        tabix -f src/\${f}.QC.mask.filt.gubbins.snps.vcf.gz

        # Concensus
        bcftools consensus \\
            -f ${params.snippy_reference_fa} \\
            src/\${f}.QC.mask.filt.gubbins.snps.vcf.gz \\
            > src/\${f}.QC.mask.filt.gubbins.snps.lineage.consensus.fa

        sed -i "s@>NC_000962.3@>\${sampleid}@g" \\
           src/\${f}.QC.mask.filt.gubbins.snps.consensus.fa

        # Copy over new vcfs and consensus files to sample directory
        cp src/\${f}.QC.mask.filt.gubbins.lineage.vcf.gz \${src_dir}/
        cp src/\${f}.QC.mask.filt.gubbins.snps.lineage.consensus.fa \${src_dir}/

    done < samples.list

    ###############################################################
    # Concatenate the consensus genomes
    ###############################################################
    cat src/*.QC.mask.filt.gubbins.snps.lineage-consensus.fa \\
        > ${lineage}.QC.mask.filt.gubbins.snps.consensus.fa

    gzip --best ${lineage}.QC.mask.filt.gubbins.snps.consensus.fa

    ###############################################################
    # Calculate the pairwise distances between the consensus genomes
    ###############################################################
    pdsm --threads ${task.cpus} --sort --progress \\
        ${lineage}.QC.mask.filt.gubbins.snps.consensus.fa.gz \\
        > ${lineage}.QC.mask.filt.gubbins.snps.consensus.pdsm.csv

    pdsm --threads ${task.cpus} --sort --progress \\
        gubbins-input.fa \\
        > ${lineage}.QC.mask.filt.snps.consensus.pdsm.csv
    """

}

/*
@author: Poppy J Hesketh Best
@date: 2026-05-06
@version: 1.0.1
@function:
    This process performs the core SNP analysis for multiple samples
        using the Snippy-core tool from the Snippy pipeline.
@details:
@references: 
    https://bitsandbugs.org/2019/11/06/two-easy-ways-to-run-iq-tree-with-the-correct-number-of-constant-sites/
@changelog
    v1.0.0-2026-05-06: Initial version
    v1.0.1-2026-05-08: Added masking with the ${param.snippy_masking}zz
*/