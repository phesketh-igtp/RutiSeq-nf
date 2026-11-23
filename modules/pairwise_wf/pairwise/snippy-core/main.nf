process SNIPPY_CORE {

    conda params.gubbins_env
    
    publishDir "${params.outDir}/db/comparison/snippy/core/", 
        mode: 'copy',
        overwrite: true

    input:
        val(sampleID_list)  // Collection of VCF files from snippy runs
        path(analysis_summary) // purely to trigger execution when pairwise analysis is run

    output:
        path("gubbins/")
        path("MTB.QC.mask.filt.snps.consensus.fa.gz")
        path("MTB.QC.mask.filt.snps.consensus.pdsm.csv")
        path("MTB.QC.mask.filt.gubbins.snps.consensus.fa.gz"), emit: core_phylo_ch
        path("MTB.QC.mask.filt.gubbins.snps.consensus.pdsm.csv")

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
        > MTB.QC.mask.filt.snps.consensus.fa.gz

    gunzip --keep MTB.QC.mask.filt.snps.consensus.fa.gz

    ###############################################################
    # Run gubbins
    ###############################################################
    run_gubbins.py \\
        --extensive-search \\
        --prefix MTB.gubbins \\
        --threads ${task.cpus} \\
        MTB.QC.mask.filt.snps.consensus.fa

    gff2bed < MTB.gubbins.recombination_predictions.gff \\
        > MTB.gubbins.recombination_predictions.bed

    sed -i 's@SEQUENCE@NC_000962.3@g' \\
        MTB.gubbins.recombination_predictions.bed
    
    # Clean up
    mkdir gubbins/; mv MTB.gubbins.* gubbins/

    ###############################################################
    # Filter VCFs with gubbins masking sites and regen consensus
    ###############################################################
    # Mask with gubbins
    for vcf in *.QC.mask.filt.vcf.gz; do

        sampleid=\$(basename ${vcf} .QC.mask.filt.vcf.gz)

        bcftools filter \\
            \${vcf} \\
            --mask-file gubbins/MTB.gubbins.recombination_predictions.gff \\
            --soft-filter gubbins \
            -Oz \\
            -o "\${sampleid}.QC.mask.filt.gubbins.vcf.gz"

        # Apply only SNPs to consensus genome
        bcftools view \\
            -f PASS \\
            -v snps \\
            "\${sampleid}.QC.mask.filt.gubbins.vcf.gz" \\
            -Oz \\
            -o "\${sampleid}.QC.mask.filt.gubbins.snps.vcf.gz"
        tabix -f "\${sampleid}.QC.mask.filt.gubbins.snps.vcf.gz"

        # Concensus
        bcftools consensus \
            -f ${params.snippy_reference_fa} \\
            "\${sampleid}.QC.mask.filt.gubbins.snps.vcf.gz" \\
            > "\${sampleid}.QC.mask.filt.gubbins.snps.consensus.fa"

        sed -i "s@>NC_000962.3@>\${sampleid}@g" \\
           "\${sampleid}.QC.mask.filt.gubbins.snps.consensus.fa"

    done

    ###############################################################
    # Concatenate the consensus genomes
    ###############################################################
    cat *.QC.mask.filt.gubbins.snps.consensus.fa \\
        > MTB.QC.mask.filt.gubbins.snps.consensus.fa
    gzip --best MTB.QC.mask.filt.gubbins.snps.consensus.fa

    ###############################################################
    # Calculate the pairwise distances between the consensus genomes
    ###############################################################
    pdsm --threads ${task.cpus} --sort --progress \\
        MTB.QC.mask.filt.gubbins.snps.consensus.fa.gz \\
        > MTB.QC.mask.filt.gubbins.snps.consensus.pdsm.csv

    pdsm --threads ${task.cpus} --sort --progress \\
        gubbins-input.fa \\
        > MTB.QC.mask.filt.snps.consensus.pdsm.csv
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-11-01
@version: 1.0.1
@function:
    This process performs the core SNP analysis for multiple samples
        using the Snippy-core tool from the Snippy pipeline.
@details:
@references: 
    https://bitsandbugs.org/2019/11/06/two-easy-ways-to-run-iq-tree-with-the-correct-number-of-constant-sites/
@changelog
    v1.0.0-2025-11-01: Initial version
    v1.0.1-2025-11-27: Updated to output variant and invariant site alignments for phylogeny
*/