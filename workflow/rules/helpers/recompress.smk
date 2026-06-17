rule helpers__recompress:
    """Base bgzip recompression rule"""
    input:
        gz="mock.fa.gz",
    output:
        gz="recompress.fa.gz",
    log:
        "recompress.log",
    conda:
        "../../environments/reference.yml"
    cache: "omit-software"
    threads: 8
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.gz} \
        | bgzip \
            --threads {threads} \
        > {output.gz} \
        ) 2> {log}
        """
