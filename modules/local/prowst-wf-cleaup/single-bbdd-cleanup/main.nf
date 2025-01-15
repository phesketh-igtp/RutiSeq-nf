process POST_SINGLE_BBDD_CLEANUP {

    input:
        tuple val(sampleID)

    script:

        """
            # rm the files that were carried over by the tuples

            # MTBseq outputs
            rm ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/*R1*fastq.gz \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/*R2*fastq.gz \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/tbdb-${sampleID}.results.txt \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/who-${sampleID}.results.txt \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mapping_and_Variant_Statistics.tab \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Strain_Classification.tab

            # SNP-Profile outputs
            rm ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/*R1*fastq.gz \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/*R2*fastq.gz \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/tbdb-${sampleID}.results.txt \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/who-${sampleID}.results.txt \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/Mapping_and_Variant_Statistics.tab \\
                ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/Strain_Classification.tab

            # TB-Profiler outputs
            rm ${params.outdir}/bbdd/tbprofiler/*R1*fastq.gz \\
                ${params.outdir}/bbdd/tbprofiler/*R2*fastq.gz \\
                ${params.outdir}/bbdd/tbprofiler/tbdb-${sampleID}.results.txt \\
                ${params.outdir}/bbdd/tbprofiler/who-${sampleID}.results.txt \\
                ${params.outdir}/bbdd/tbprofiler/Mapping_and_Variant_Statistics.tab \\
                ${params.outdir}/bbdd/tbprofiler/Strain_Classification.tab

            # Readsqc outputs - remove mtc reads and any lingering reads
            rm ${params.outdir}/bbdd/read-qc/mtbc_reads/mtbc-${sampleID}_R1.fastq.gz \\
                ${params.outdir}/bbdd/read-qc/mtbc_reads/mtbc-${sampleID}_R2.fastq.gz \\
                ${params.outdir}/bbdd/read-qc/mtbc_reads/*R1*fastq.gz \\
                ${params.outdir}/bbdd/read-qc/mtbc_reads/*R2*fastq.gz

            # compress the outputs from MTBSeq mpileup
            gzip --best ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.mpileup
            gzip --best ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.mpileuplog

        """

}