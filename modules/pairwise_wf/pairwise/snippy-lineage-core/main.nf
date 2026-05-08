process SNIPPY_LINEAGE_CORE {

    tag "${lineage}"

    conda params.snippy_env
    
    publishDir "${params.outDir}/db/comparison/snippy/${lineage}/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple val(lineage), 
            val(sampleIDs), // Collection of VCF files from snippy runs
            val(sampleID_count)

    output:
        // unmasked outputs
        path("${lineage}-core.aln")
        path("${lineage}-core.tab")
        path("${lineage}-core.full.aln")
        path("${lineage}-core.txt")
        path("${lineage}-core.vcf")
        path("${lineage}-core.mat.tsv")
        // masked oututs
        path("${lineage}-masked.aln")
        path("${lineage}-masked.tab")
        path("${lineage}-masked.full.aln")
        path("${lineage}-masked.txt")
        path("${lineage}-masked.vcf")
        path("${lineage}-masked.mat.tsv")
        // Reference
        path("${lineage}.ref.fa")

    script:

    """
    # create the list of the sampleIDs within that lineage
            echo "${sampleIDs.join('\n')}" | sort | uniq > samples.list

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
    snippy_directories=\$(echo snippyDir_*)

    # The main alignments/core
        snippy-core \\
            --prefix ${lineage}-core \\
            --ref ${params.snippy_reference} \\
            \$snippy_directories

        seqkit seq \\
            --line-width 0 \\
            --threads ${task.cpus} \\
            ${lineage}-core.full.aln \\
            > tmp.core.full.aln
        mv tmp.core.full.aln ${lineage}-core.full.aln

    # The main masked alignments/core
        snippy-core \\
            --prefix ${lineage}-masked \\
            --ref ${params.snippy_reference} \\
            --mask ${params.snippy_masking} \\
            \$snippy_directories

        seqkit seq \\
            --line-width 0 \\
            --threads ${task.cpus} \\
            ${lineage}-masked.full.aln \\
            > tmp.${lineage}-core.m
        mv tmp.${lineage}-core.m ${lineage}-masked.full.aln

    # Get snp distances
    snp-dists \\
        -k -j ${task.cpus} \\
        ${lineage}-core.aln \\
        > ${lineage}-core.mat.tsv
    
    snp-dists \\
        -k -j ${task.cpus} \\
        ${lineage}-masked.aln \\
        > ${lineage}-masked.mat.tsv

    # Clean up
        # Remove the prefix
        sed -i 's@snippyDir_@@g' ${lineage}-core*
        mv ${lineage}-core.ref.fa > ${lineage}.ref.fa
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
    v1.0.1-2026-05-08: Added masking with the ${param.snippy_masking}
*/