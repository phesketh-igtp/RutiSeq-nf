process TBPROFILER_PROFILE_WHO {


    conda "bioconda::mtbseq=1.1.0"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mtbseq:1.1.0--hdfd78af_0' :
        'biocontainers/mtbseq:1.1.0--hdfd78af_0' }"

    input:

    output:

    script:

        """

        """

    stub:
        """

        """

}