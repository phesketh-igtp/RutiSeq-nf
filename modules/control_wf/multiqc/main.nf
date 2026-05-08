process MULTIQC {

    conda params.readQC_env

    storeDir "${params.outDir}/db/qc/${params.runID}/"

    input:
        path(samplesheet, stageAs: 'samplesheet.csv')

    output:
        path("multiQC-out/multiqc_report.html")

    script:
    """
    mkdir -p reads/ multiQC-out/ fastqc/
    
    # Create symbolic links to the fastq files
    while IFS=',' read -r identities samples fastq_1 fastq_2 type; do

        fastqc \\
            \${fastq_1} \\
            \${fastq_2} \\
            --threads ${task.cpus} \\
            --outdir fastqc/

    done < <(tail -n +2 samplesheet.csv)

    # MultiQC report
    multiqc fastqc/* --outdir multiQC-out/
    """
}


/*
@author: Poppy J hesketh Best
@date: 2025-11-17
@version: v1.0.0
@description:
    Runs fastQC then multiQC on the reads. Saves only the multiQC HTML
@changelog:
    v1.0.0-2025-11-17: Functioning module created.
*/
