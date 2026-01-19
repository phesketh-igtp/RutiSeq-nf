process VERSION_LOGGING {

    publishDir "${params.outDir}/results/${params.runID}/", 
        mode: 'copy',
        overwrite: true

    input:
        val(runID)

    output:
        file("${params.runID}_version_log.txt"), emit: version_logs

    script:
    """
    # concatenate yml files from conda
    for yaml in $params.envs_dir/*.yml; do

        sed '1,/^dependencies:/d' \$yaml > \$(basename \$yaml).version_log.txt

    done

        for yaml in $params.envs_dir/*.yaml; do

        sed '1,/^dependencies:/\n/d' \$yaml > \$(basename \$yaml).version_log.txt

    done

    # reformat into a txt file
    cat *.version_log.txt \\
        | sort version_log.yml \\
        | uniq \\
        | sed 's@  - @@g' \\
        | sed 's@=@\t@g' \\
        | sed 's@::@\t@g'> ${params.runID}_version_log.txt
    """
}