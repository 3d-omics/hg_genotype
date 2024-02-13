rule index:
    """Build bowtie2 index

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=REFERENCE / "genome.fa.gz",
    output:
        multiext(
            f"{INDEX}/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        INDEX / "build.log",
    conda:
        "__environment__.yml"
    params:
        output_path=INDEX / "genome",
        extra=params["bowtie2"]["extra"],
    threads: 8
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.reference} \
            {params.output_path} \
        2> {log} 1>&2
        """
