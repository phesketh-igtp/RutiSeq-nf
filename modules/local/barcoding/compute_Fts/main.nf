process COMPUTE_FTS {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-22
    @version: 1.0.0
    @description: 
        This process calcualte the Fst value from the cluster VCF using the
        merged lineage VCFs
    @changelog:
        v1.0.0-2025-04-22:  Initial version
*/

    tag "${clusterID}:${lineage}"

    publishDir "${params.outDir}/barcoding/tabs/", mode: 'copy', overwrite: true

    conda params.snp_profiling_env
    
    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_snp_profiling
            else if (workflow.containerEngine == 'docker') return params.docker_snp_profiling
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_snp_profiling
            else return null
        }

    input:
        tuple val(clusterID),
            val(lineage),
            val(sampleIDs), 
            path(lineagePopList), 
            path(mergedVcf), 
            path(mergedVcfIndex)


    output:
        tuple val(clusterID),
            val(lineage),
            path("${clusterID}.merge.vcf.gz"),
            path("${clusterID}.merge.vcf.tab")

    script:
    
        """
        # create a list of the sampleIDs of the clsuterID
            echo "${sampleIDs.join('\n')}" > ${clusterID}-pop.list

        # Unzip the VCF as vcftools doesnt work on zipped ones
            bgzip -d ${mergedVcf} > merged.vcf

        # calculate Fst statistics between individuals of different populations
            vcftools --vcf merged.vcf \\
                    --weir-fst-pop ${lineagePopList} \\
                    --weir-fst-pop ${clusterID}.pop.list \\
                    --out ${clusterID}_vs_${lineage}.tsv
        
        # rm the unzipped vcf file
            rm merged.vcf
        """

}