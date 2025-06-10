process POST_SINGLE_BBDD_CLEANUP {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process is responsible for cleaning up the output files generated
        during the analysis of the pairwise clusters. It removes unnecessary files
        and compresses the mpileup files to save space.
*/

    tag "${runID}"

    input:
        val runID
        path pairwise_clusters 

    script:

        """
        # Check and remove files only if they exist

        for file in \\
            "${params.outDir}/bbdd/mtbseq/samples/*/*R1*fastq.gz" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/*R2*fastq.gz" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/tbdb-*.results.txt" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/who-*.results.txt" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/Mapping_and_Variant_Statistics.tab" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/Strain_Classification.tab" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/SNP-Profiles/*R1*fastq.gz" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/SNP-Profiles/*R2*fastq.gz" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/SNP-Profiles/tbdb-*.results.txt" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/SNP-Profiles/who-*.results.txt" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/SNP-Profiles/Mapping_and_Variant_Statistics.tab" \\
            "${params.outDir}/bbdd/mtbseq/samples/*/SNP-Profiles/Strain_Classification.tab" \\
            "${params.outDir}/bbdd/tbprofiler/*_mtbc_R1.fastq.gz" \\
            "${params.outDir}/bbdd/tbprofiler/*_mtbc_R2.fastq.gz" \\
            "${params.outDir}/bbdd/tbprofiler/tbdb-*.results.txt" \\
            "${params.outDir}/bbdd/tbprofiler/who-only/*_mtbc_R1.fastq.gz" \\
            "${params.outDir}/bbdd/tbprofiler/who-only/*_mtbc_R2.fastq.gz" \\
            "${params.outDir}/bbdd/tbprofiler/who-only/who-*.results.txt" \\
            "${params.outDir}/bbdd/tbprofiler/who-only/tbdb-*.results.txt" \\
            "${params.outDir}/bbdd/tbprofiler/Mapping_and_Variant_Statistics.tab" \\
            "${params.outDir}/bbdd/tbprofiler/Strain_Classification.tab" \\
            "${params.outDir}/bbdd/read-qc/mtbc_reads/*_mtbc_R1.fastq.gz" \\
            "${params.outDir}/bbdd/read-qc/mtbc_reads/*_mtbc_R2.fastq.gz"; 
        do
            if [ -f \$file ] || [ -e \$file ]; then rm \$file; fi
        done

        # Compress the outputs from MTBSeq mpileup
        if [ -f "${params.outDir}/bbdd/mtbseq/samples/*/Mpileup/*.mpileup" ]; then
            gzip --best "${params.outDir}/bbdd/mtbseq/samples/*/Mpileup/*.mpileup"
        fi

        if [ -f "${params.outDir}/bbdd/mtbseq/samples/*/Mpileup/*.mpileuplog" ]; then
            gzip --best "${params.outDir}/bbdd/mtbseq/samples/*/Mpileup/*.mpileuplog"
        fi

        """

}