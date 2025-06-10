process POST_SUMMARY_CLEANUP {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process performs cleanup of the handover files and compresses the outputs from MTBSeq mpileup.
        It removes any existing files in the specified locations and compresses the FASTA files in the output directory.
*/

    input:
        path(handover)

    script:

        """
        # Check and remove files only if they exist

        for file in \\
            "${params.outDir}/results/snps/cleanup-handover" \\
            "${params.outDir}/results/networks/*_joint_*phylo.tab"
        do
            if [ -f \$file ] || [ -e \$file ]; then rm \$file; fi
        done

        # Compress the outputs from MTBSeq mpileup
        for fasta in "${params.outDir}/bbdd/results/networks/fasta/*"; do
            if [ -f \${fasta} ]; then gzip --best \${fasta}; fi
        """

}