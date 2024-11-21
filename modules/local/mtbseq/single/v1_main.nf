process MTBSEQ_SINGLE {

    conda "bioconda::mtbseq=1.1.0"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mtbseq:1.1.0--hdfd78af_0' :
        'biocontainers/mtbseq:1.1.0--hdfd78af_0' }"

    input:
        tuple val(old_name), val(sampleID), path(forward), path(reverse)
        path "${sampleID}"
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        path "lineage-compile/Amend", emit: called // paths defined by the sampleID
        path "lineage-compile/Position_Tables", emit: position_tables_dir
        path "lineage-compile/Classification", emit: classification_dir
        path "lineage-compile/Statistics", emit: statistics_dir
        path "lineage-compile/Statistics/Mapping_and_Variant_Statistics.tab", emit: statistics
        path "lineage-compile/Classification/Strain_Classification.tab", emit: classification
        path "lineage-compile/Called/*gatk_position_variants*.tab", emit: position_variants
        path "lineage-compile/Position_Tables/*.gatk_position_table.tab", emit: position_tables
        path "versions.yml", emit: versions

    script:
        def args = task.ext.args ?: " --minbqual ${params.minbqual} --mincovf ${params.mincovf} --mincovr ${params.mincovr} --minphred ${params.minphred} --minfreq ${params.minfreq} --unambig ${params.unambig} --window ${params.window}"

        """

        #·············  PREPARE WORING DIRECTORY FOR PROCESSING  ·············#

            mkdir -p lineage-compile
 
            # Create symbolic link to data in individual directory
            ln -s ${forward} ${sampleID}_R1.fastq.gz
            ln -s ${reverse} ${sampleID}_R2.fastq.gz

            # Generate MTBSeq sample.txt (this might be redundant)
            ls ${forward} | sed 's/.*\///' | sed 's/_R1\.fastq\.gz//g' | sed 's/_/\t/' | cut -f1 > tmp.sampleID.col1
            ls ${forward} | sed 's/.*\///' | sed 's/_R1\.fastq\.gz//g' | sed 's/_/\t/' | cut -f2  > tmp.libraryID.col2
            paste tmp.sampleID.col1 tmp.libraryID.col2 > ${sampleID}_sample.txt; rm tmp.sampleID.col1 tmp.libraryID.col2

        #·············  MTBSEQ  ·············#

            MTBseq --step TBfull \
                --samples ${sampleID}_sample.txt
                --project ${sampleID}
                --thread ${task.cpus} \
                ${args}

        """

    stub:
        """
        sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

        echo "MTBseq --step TBfull \
            --thread ${task.cpus} \
            --project ${params.project} \
            --minbqual ${params.minbqual} \
            --mincovf ${params.mincovf} \
            --mincovr ${params.mincovr} \
            --minphred ${params.minphred} \
            --minfreq ${params.minfreq} \
            --unambig ${params.unambig} \
            --window ${params.window} \
            --distance ${params.distance} \
            --resilist ${ref_resistance_list} \
            --intregions ${ref_interesting_regions} \
            --categories ${ref_gene_categories} \
            --basecalib ${ref_base_quality_recalibration} "

        mkdir GATK_Bam
        touch GATK_Bam/stub.gatk.bam
        touch GATK_Bam/stub.gatk.bai
        touch GATK_Bam/stub.gatk.bamlog
        touch GATK_Bam/stub.gatk.grp
        touch GATK_Bam/stub.gatk.intervals
        mkdir Bam
        touch Bam/stub.bam
        touch Bam/stub.bai
        touch Bam/stub.bamlog
        mkdir Called
        touch Called/stub.gatk_position_uncovered_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
        touch Called/stub.gatk_position_variants_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
        mkdir Mpileup
        touch Mpileup/stub.gatk.mpileup
        touch Mpileup/stub.gatk.mpileuplog
        mkdir Classification
        touch Classification/Strain_Classification.tab
        mkdir Position_Tables
        touch Position_Tables/stub.gatk_position_table.tab
        mkdir Statistics
        touch Statistics/Mapping_and_Variant_Statistics.tab

        """

}