rule align__index__bwa:
    """Build genome index with bwa"""
    input:
        reference=REFERENCE / "genome.fa.gz",
    output:
        multiext(f"{INDEX}/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
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
        bwa index \
            -p {params.output_path} \
            {input.reference} \
            {params.extra} \
        2> {log} 1>&2
        """
