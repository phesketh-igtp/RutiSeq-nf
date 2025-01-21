process POST_SINGLE_BBDD_CLEANUP {

    tag "${runID}"

    input:
        val runID
        path pairwise_clusters 

    script:

        """
        # Check and remove files only if they exist

        for file in \\
            "${params.outdir}/bbdd/mtbseq/samples/*/*R1*fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/*R2*fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/tbdb-*.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/who-*.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/Mapping_and_Variant_Statistics.tab" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/Strain_Classification.tab" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/SNP-Profiles/*R1*fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/SNP-Profiles/*R2*fastq.gz" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/SNP-Profiles/tbdb-*.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/SNP-Profiles/who-*.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/SNP-Profiles/Mapping_and_Variant_Statistics.tab" \\
            "${params.outdir}/bbdd/mtbseq/samples/*/SNP-Profiles/Strain_Classification.tab" \\
            "${params.outdir}/bbdd/tbprofiler/*_mtbc_R1.fastq.gz" \\
            "${params.outdir}/bbdd/tbprofiler/*_mtbc_R2.fastq.gz" \\
            "${params.outdir}/bbdd/tbprofiler/tbdb-*.results.txt" \\
            "${params.outdir}/bbdd/tbprofiler/who-only/*_mtbc_R1.fastq.gz" \\
            "${params.outdir}/bbdd/tbprofiler/who-only/*_mtbc_R2.fastq.gz" \\
            "${params.outdir}/bbdd/tbprofiler/who-only/who-*.results.txt" \\
            "${params.outdir}/bbdd/tbprofiler/who-only/tbdb-*.results.txt" \\
            "${params.outdir}/bbdd/tbprofiler/Mapping_and_Variant_Statistics.tab" \\
            "${params.outdir}/bbdd/tbprofiler/Strain_Classification.tab" \\
            "${params.outdir}/bbdd/read-qc/mtbc_reads/*_mtbc_R1.fastq.gz" \\
            "${params.outdir}/bbdd/read-qc/mtbc_reads/*_mtbc_R2.fastq.gz"; 
        do
            if [ -f \$file ] || [ -e \$file ]; then rm \$file; fi
        done

        # Compress the outputs from MTBSeq mpileup
        if [ -f "${params.outdir}/bbdd/mtbseq/samples/*/Mpileup/*.mpileup" ]; then
            gzip --best "${params.outdir}/bbdd/mtbseq/samples/*/Mpileup/*.mpileup"
        fi

        if [ -f "${params.outdir}/bbdd/mtbseq/samples/*/Mpileup/*.mpileuplog" ]; then
            gzip --best "${params.outdir}/bbdd/mtbseq/samples/*/Mpileup/*.mpileuplog"
        fi

        """

}