process ENSURE_COMPLETION_AND_HANDOFF {
    input:
    val completed_processes

    output:
    path "single_workflow_complete.txt", emit: handoff_signal

    script:
    """
    echo "Single workflow completed at \$(date)" > single_workflow_complete.txt
    echo "Completed processes: ${completed_processes}" >> single_workflow_complete.txt
    """
}