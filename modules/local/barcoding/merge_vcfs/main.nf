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
            path(path_to_vcfs),

    output:
        tuple val(lineage),
            path("${lineage}.merged.vcf.gz"),
            path("${lineage}.merged.vcf.gz.tbi"),

    script:
    
        """
        # Use a file to access all the VCF files of genomes within a single lineage
            bcftools merge \\
                –file-list ${path_to_vcfs} \\
                -Oz -o ${lineage}.merged.vcf.gz

        # create the index file
            bcftools index -t ${lineage}.merged.vcf.gz
        """

}