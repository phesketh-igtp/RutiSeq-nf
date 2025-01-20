process POST_SINGLE_BBDD_CLEANUP {

    tag "${sampleID}"

    input:
        tuple val(sampleID)

    script:

        """
        # Check and remove files only if they exist

        for file in \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/*R1*fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/*R2*fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/tbdb-${sampleID}.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/who-${sampleID}.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mapping_and_Variant_Statistics.tab" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Strain_Classification.tab" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/*R1*fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/*R2*fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/tbdb-${sampleID}.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/who-${sampleID}.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/Mapping_and_Variant_Statistics.tab" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/Strain_Classification.tab" \\
            "${params.outdir}/bbdd/tbprofiler/${sampleID}_mtbc_R1.fastq.gz" \\
            "${params.outdir}/bbdd/tbprofiler/${sampleID}_mtbc_R2.fastq.gz" \\
            "${params.outdir}/bbdd/tbprofiler/tbdb-${sampleID}.results.txt" \\
            "${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}_mtbc_R1.fastq.gz" \\
            "${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}_mtbc_R2.fastq.gz" \\
            "${params.outdir}/bbdd/tbprofiler/who-only/who-${sampleID}.results.txt" \\
            "${params.outdir}/bbdd/tbprofiler/Mapping_and_Variant_Statistics.tab" \\
            "${params.outdir}/bbdd/tbprofiler/Strain_Classification.tab" \\
            "${params.outdir}/bbdd/read-qc/mtbc_reads/${sampleID}_mtbc_R1.fastq.gz" \\
            "${params.outdir}/bbdd/read-qc/mtbc_reads/${sampleID}_mtbc_R2.fastq.gz"; 
        do
            if [ -f \$file ] || [ -e \$file ]; then rm \$file; fi
        done

        # Compress the outputs from MTBSeq mpileup
        if [ -f "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.mpileup" ]; then
            gzip --best "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.mpileup"
        fi

        if [ -f "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.mpileuplog" ]; then
            gzip --best "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.mpileuplog"
        fi

        """

}