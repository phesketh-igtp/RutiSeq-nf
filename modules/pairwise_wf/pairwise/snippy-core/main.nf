process SNIPPY_CORE {

    conda params.snippy_env
    
    publishDir "${params.outDir}/db/comparison/snippy/", mode: 'copy'

    input:
        val(sampleID_list)  // Collection of VCF files from snippy runs
        path(pairwise_analysis_list) // purely to trigger execution when pairwise analysis is run

    output:
        path("core.snps")
        path("core.vcf")
        path("core.aln")
        path("core.aln.full")
        path("core.tab")
        path("core.txt")
        path("core.masked.distance.mat")
        path("core.masked.distance.mat.tsv")

        path("core.masked.snps"),       emit: snippy_core_snps
        path("core.masked.vcf"),        emit: snippy_core_vcf
        path("core.masked.aln"),        emit: snippy_core_alignment
        path("core.masked.aln.full"),   emit: snippy_core_alignment_full
        path("core.masked.distance.mat")
        path("core.masked.distance.mat.tsv")
        path("core.masked.tab")
        path("core.masked.txt")

        tuple path("phylo.variant-sites.aln"),
            path("phylo.invariant-sites.aln"), emit: snippy_core_phylo_alignment

        path("samples.list")

    script:

    """
    # Create the sample list of available results
    ls ${params.outDir}/db/samples/ > samples.list

    # Assign the base path
    BASE="${params.outDir}/db/samples"

    while read -r f; do
        src_dir="\$BASE/\$f/snippy"
        dest_dir="snippyDir_\$f"

        # Check source directory
        if [[ ! -d "\$src_dir" ]]; then
            echo "No snippy directory for sample: \$f"
            continue
        fi

        # Check for target files
        vcf_files=("\$src_dir"/*.vcf)
        fa_files=("\$src_dir"/*.aligned.fa)

        vcf_exists=false
        fa_exists=false

        [[ -e "\${vcf_files[0]}" ]] && vcf_exists=true
        [[ -e "\${fa_files[0]}" ]] && fa_exists=true

        # Skip sample entirely if neither file type exists
        if ! \$vcf_exists && ! \$fa_exists; then
            echo "No VCF or aligned FASTA files found for \$f — skipping."
            continue
        fi

        # Only create dest directory if needed
        mkdir -p "\$dest_dir"

        # Link only existing files
        if \$vcf_exists; then
            ln -sf "\$src_dir"/*.vcf "\$dest_dir"/
            echo "Linked VCF files for \$f"
        fi

        if \$fa_exists; then
            ln -sf "\$src_dir"/*.aligned.fa "\$dest_dir"/
            echo "Linked aligned FASTA files for \$f"
        fi

    done < samples.list

    # Create path of directories for snippy core
    PATH=\$(echo snippyDir_*)

    # The main alignments/core
        snippy-core --prefix core \\
            --ref ${params.snippy_reference} \\
            \$PATH

        seqkit seq -w 0 core.full.aln > tmp.core.full.aln
        mv tmp.core.full.aln core.full.aln
        snp-dists -j 8 core.full.aln > core.distance.mat
        snp-dists -j 8 -m core.full.aln > core.distance.mat.tsv

    # The main masked alignments/core
        snippy-core --prefix core.masked \\
            --ref ${params.snippy_reference} \\
            --mask ${params.snippy_reference} \\
            \$PATH

        seqkit seq -w 0 core.masked.full.aln > tmp.core.masked.full.aln
        mv tmp.core.masked.full.aln core.masked.full.aln
        snp-dists -j 8 core.masked.full.aln > core.masked.distance.mat
        snp-dists -j 8 -m core.masked.full.aln > core.masked.distance.mat.tsv

    # Remove the prefix
        sed -i 's@snippyDir_@@g' core.*
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