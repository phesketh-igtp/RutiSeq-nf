process POST_SINGLE_BBDD_CLEANUP {

    /*
        This module clean up all the intermediate files that are carried over through the process
            of emits/publications - this has been an awkward workaround to ensuring that one 
            genome that do not have the necessary outputs for the PAIRWISE_WF() are analysed.
        Unfortunately at each module step then the tuple containing all the paths to the previous
            are published in the output directory (mode: 'copy') - this final module just ensures that 
            all intermediate files are removed from the publish directory - even though they would have
            been remove during the process.
    */

    tag "${sampleID}"

    array 100

    input:
        tuple val(sampleID)

    script:

        """
        # Check and remove files only if they exist

        for file in \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/${sampleID}_mtbc_R1.fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/${sampleID}_mtbc_R2.fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/tbdb-${sampleID}.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/who-${sampleID}.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/Mapping_and_Variant_Statistics.tab" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/Strain_Classification.tab" \\
        do
            if [ -f "\${file}" ] || [ -e "\${file}" ]; then rm "\${file}"; fi
        done

        # Compress the outputs from MTBSeq mpileup
            if [ -f "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.gatk.mpileup" ]; then
                gzip --force --best "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.gatk.mpileup"
            fi

            if [ -f "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.gatk.mpileuplog" ]; then
                gzip --force --best "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.gatk.mpileuplog"
            fi

        # Compress the outputs from MTBSeq mpileup
            if [ -f "${params.outdir}/bbdd/read-qc/tables/${sampleID}.kaiju.out" ]; then
                gzip --force --best "${params.outdir}/bbdd/read-qc/tables/${sampleID}.kaiju.out"
            fi

        """

}