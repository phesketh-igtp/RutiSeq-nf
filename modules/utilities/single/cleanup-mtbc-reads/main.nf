process CLEANUP_MTBC_READS {
    input:
    tuple path(r1), path(r2)

    script:
    """
    rm -f $r1 $r2
    """
}