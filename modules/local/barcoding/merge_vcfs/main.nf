process MERGE_VCFS {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-22
    @version: 1.0.0
    @description: 
        This process merges the VCFs from all the sample in the same lineage
    @changelog:
        v1.0.0-2025-04-22:  Initial version
*/

    tag "${lineage}"


    conda params.snp_profiling_env
    
    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_snp_profiling
            else if (workflow.containerEngine == 'docker') return params.docker_snp_profiling
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_snp_profiling
            else return null
        }

    input:
        tuple val(lineage), 
            val(clusterID), 
            val(sampleIDs), 
            path(vcf_paths)

    output:
        tuple val(clusterID),
            val(lineage),
            val(sampleIDs),
            path("${lineage}.pop.list"),
            path("${lineage}.merged.vcf.gz"),
            path("${lineage}.merged.vcf.gz.tbi"), emit: merged_vcs_tuple

    script:
    
        """
        # Use a file to access all the VCF files of genomes within a single lineage
            bcftools merge \\
                -file-list ${vcf_paths} \\
                -Oz -o ${lineage}.merged.vcf.gz

        # create the index file
            bcftools index -t ${lineage}.merged.vcf.gz

        # Create a list of the sampleID within the lineage
            echo "${sampleIDs.join('\n')}" > ${lineage}.pop.list
        """

}