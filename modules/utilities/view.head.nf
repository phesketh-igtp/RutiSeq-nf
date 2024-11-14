process VIEW_HEAD {
    input:
    path file

    output:
    stdout

    script:
    """
    echo "First 10 lines of ${file}:"
    head -n 10 ${file}
    """
}