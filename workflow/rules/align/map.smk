rule align__map__bwamem__:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=READS / "{sample}.{library}_1.fq.gz",
        reverse_=READS / "{sample}.{library}_2.fq.gz",
        idx=multiext(f"{INDEX}/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        reference=REFERENCE / "genome.fa.gz",
    output:
        cram=MAP / "{sample}.{library}.cram",
    log:
        MAP / "{sample}.{library}.log",
    params:
        index_prefix=INDEX / "genome",
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        read_group_header=compose_read_group_header,
    threads: 24
    conda:
        "__environment__.yml"
    resources:
        mem_mb=30000,
        runtime=1440,
    shell:
        """
        (bwa mem \
            -t {threads} \
            -R '{params.read_group_header}' \
            {params.index_prefix} \
            {params.extra} \
            {input.forward_} \
            {input.reverse_} \
        | samtools sort \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule align__map__bwamem__all:
    """Collect the results of `bowtie2_map_one` for all libraries"""
    input:
        [MAP / f"{sample}.{library}.cram" for sample, library in SAMPLE_LIB],
