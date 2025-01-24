process POST_SUMMARY_CLEANUP {

    input:
        path(handover)

    script:

        """
        # Check and remove files only if they exist

        for file in \\
            "${params.outdir}/results/snps/cleanup-handover" \\
            "${params.outdir}/results/networks/*_joint_*phylo.tab"
        do
            if [ -f \$file ] || [ -e \$file ]; then rm \$file; fi
        done

        # Compress the outputs from MTBSeq mpileup
        for fasta in "${params.outdir}/bbdd/results/networks/fasta/*"; do
            if [ -f \${fasta} ]; then gzip --best \${fasta}; fi
        """

}