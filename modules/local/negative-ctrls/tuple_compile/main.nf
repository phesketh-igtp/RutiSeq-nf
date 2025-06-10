process CN_COMPILE_SUMMARY_TUPLE {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-07
    @version: 1.0.0
    @description:
        This process takes the output from the Kraken2, Tb-Profiler and MTBseq processes
        and compiles them into a single tuple for the negative control analysis.
    @changelog:
        2025-04-07: Created process.
*/

    tag "$sampleID"

    array 100

    input:
    tuple val(sampleID), path(files)

    output:
    tuple val(sampleID), emit: combined_cn_tuple

    script:
    """
    """
}