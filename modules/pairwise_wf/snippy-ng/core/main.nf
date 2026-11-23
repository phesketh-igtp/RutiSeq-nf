process SNIPPY_CORE  {

    conda params.gubbins_env
    
    publishDir "${params.outDir}/db/comparison/snippy/RutiSeq/", 
        mode: 'copy',
        overwrite: true

    input:
        path(lineage_sample_tuple)

    output:
        path("gubbins/")
        path("RutiSeq.QC.mask.filt.snps.consensus.fa.gz")
        path("RutiSeq.QC.mask.filt.snps.consensus.pdsm.csv")
        path("RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.fa.gz")
        path("RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.pdsm.csv")

        path("RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.fa.gz"),
                emit: snippy_lineage_phylo_ch

    script:

    """
    ###############################################################
    # Link all the VCF files
    ###############################################################
    # Get updated list of samples in DB
    ls ${params.outDir}/db/samples/ | sort | uniq > samples.list

    # Assign the base path
    BASE="${params.outDir}/db/samples"

    while read -r f; do

        src_dir="\$BASE/\$f/snippy-ng/"
        dest_dir="src/"

        # Only create dest directory if needed
        mkdir -p src/
        ln -sf "\$src_dir"/*.QC.mask.filt.vcf.gz src/
        ln -sf "\$src_dir"/*.QC.mask.filt.snps.consensus.fa.gz src/

    done < samples.list

    ###############################################################
    # Concatenate the consensus sequences
    ###############################################################
    cat src/*.QC.mask.filt.snps.consensus.fa.gz \\
        > RutiSeq.QC.mask.filt.snps.consensus.fa.gz

    gunzip --keep RutiSeq.QC.mask.filt.snps.consensus.fa.gz

    ###############################################################
    # Run gubbins
    ###############################################################
    run_gubbins.py \\
        --extensive-search \\
        --prefix RutiSeq.gubbins \\
        --threads ${task.cpus} \\
        RutiSeq.QC.mask.filt.snps.consensus.fa

    gff2bed < RutiSeq.gubbins.recombination_predictions.gff \\
        > RutiSeq.gubbins.recombination_predictions.bed

    sed -i 's@SEQUENCE@NC_000962.3@g' \\
        RutiSeq.gubbins.recombination_predictions.bed
    
    # Clean up
    mkdir gubbins/; mv RutiSeq.gubbins.* gubbins/

    ###############################################################
    # Filter VCFs with gubbins masking sites and regen consensus
    ###############################################################
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
        cp src/\${f}.QC.mask.filt.gubbins.core.vcf.gz \${src_dir}/
        cp src/\${f}.QC.mask.filt.gubbins.snps.core.consensus.fa \${src_dir}/

    done < samples.list

    ###############################################################
    # Concatenate the consensus genomes
    ###############################################################
    cat src/*.QC.mask.filt.gubbins.snps.lineage.consensus.fa \\
        > RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.fa

    gzip --best RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.fa

    ###############################################################
    # Calculate the pairwise distances between the consensus genomes
    ###############################################################
    pdsm --threads ${task.cpus} --sort --progress \\
        RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.fa.gz \\
        > RutiSeq.QC.mask.filt.gubbins.snps.lineage.consensus.pdsm.csv

    pdsm --threads ${task.cpus} --sort --progress \\
        RutiSeq.QC.mask.filt.snps.consensus.fa \\
        > RutiSeq.QC.mask.filt.snps.consensus.pdsm.csv
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