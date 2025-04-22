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
            path(merge_vcf),
            path(lineage_pop),
            path(cluster_pop)

    output:
        tuple val(clusterID),
            val(lineage),
            path("${clusterID}.merge.vcf.gz"),
            path("${clusterID}.merge.vcf.tab")

    script:
    
        """
        # Unzip the VCF as vcftools doesnt work on zipped ones
            bgzip -d ${merge_vcf} > merged.vcf

        # calculate Fst statistics between individuals of different populations
            vcftools --vcf merged.vcf \\
                    --weir-fst-pop ${lineage_pop} \\
                    --weir-fst-pop ${cluster_pop} \\
                    --out ${clusterID}_vs_${lineage}.tsv
        
        # rm the unzipped vcf file
            rm merged.vcf
        """

}