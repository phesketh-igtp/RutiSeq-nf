process FILE_CHECK {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description: 
        This process checks for the existence of specific files related to MTBseq and TBProfiler results.
        It generates a CSV file with the sample ID and paths to the forward and reverse reads, as well 
        as the paths to the MTBseq and TBProfiler results if they exist.
*/

    tag "$sampleID"

    //array 200

    input:
        tuple val(sampleID), path(forward), path(reverse)

    output:
        path("samples.tuple.csv"), emit: sample_paths

    script:
        def forward_path = forward.toRealPath()
        def reverse_path = reverse.toRealPath()
        def mtbseq_class = "${params.outDir}/bbdd/mtbseq/samples/${sampleID}/Classification/${sampleID}.Strain_Classification.tab"
        def mtbseq_stats = "${params.outDir}/bbdd/mtbseq/samples/${sampleID}/Statistics/${sampleID}.Mapping_and_Variant_Statistics.tab"
        def mtbseq_pos   = "${params.outDir}/bbdd/mtbseq/samples/${sampleID}/Position_Tables/${sampleID}.gatk_position_table.tab"
        def mtbseq_vars  = "${params.outDir}/bbdd/mtbseq/samples/${sampleID}/Called/${sampleID}.gatk_position_variants_cf4_cr4_fr75_ph4_outmode000.tab"
        def tbdb_out     = "${params.outDir}/bbdd/tbprofiler/results/tbdb-${sampleID}.results.txt"
        def who_out      = "${params.outDir}/bbdd/tbprofiler/who-only/results/who-${sampleID}.results.txt"
        def mtbseq_vcf   = "${params.outDir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/${sampleID}.gatk.vcf.gz"

        """
        if [ -f ${mtbseq_class} ] && [ -f ${mtbseq_stats} ] && [ -f ${mtbseq_pos} ] && [ -f ${mtbseq_vars} ] && [ -f ${tbdb_out} ] && [ -f ${who_out} ] && [ -f ${mtbseq_vcf} ]; then       
            echo "${sampleID},,,${mtbseq_class},${mtbseq_stats},${mtbseq_pos},${mtbseq_vars},${tbdb_out},${who_out},${mtbseq_vcf}" > samples.tuple.csv
            echo "DEBUG: Added to samples.txt: ${sampleID}" >&2    
        else
            echo "${sampleID},${forward_path},${reverse_path},,,,,,," > samples.tuple.csv
            echo "DEBUG: Added to samples.txt: ${sampleID}" >&2
        fi

        # Ensure the process doesn't fail if one of the files doesn't exist
        touch samples.tuple.csv
        """

}