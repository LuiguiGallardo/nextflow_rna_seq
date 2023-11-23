params.str = 'your name'

process splitLetters {
    output:
    path 'chunk_*'

    """
    printf 'Hello ${params.str}' | split -b 100 - chunk_
    """
}

process convertToUpper {
    input:
    path x

    output:
    stdout

    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

workflow {
    splitLetters | flatten | convertToUpper | view { it.trim() }
}