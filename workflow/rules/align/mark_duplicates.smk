rule align__mark_duplicates__:
    """Mark duplicates in a single chromosome from a single library"""
    input:
        cram=get_crams_for_mark_duplicates,
        reference=REFERENCE / "genome.fa.gz",
    output:
        bam=pipe(MARK_DUPLICATES / "{sample}.bam"),
        metrics=MARK_DUPLICATES / "{sample}.metrics.tsv",
    log:
        MARK_DUPLICATES / "{sample}.bam.log",
    params:
        input_cram=compose_input_line_for_mark_duplicates,
    conda:
        "__environment__.yml"
    resources:
        mem_mb=8000,
        runtime=360,
    threads: 0
    shell:
        """
        gatk MarkDuplicates \
            {params.input_cram} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            --COMPRESSION_LEVEL 0 \
            --REFERENCE_SEQUENCE {input.reference} \
        2> {log} 1>&2
        """


rule align__mark_duplicates__bam_to_cram:
    input:
        bam=MARK_DUPLICATES / "{sample}.bam",
        reference=REFERENCE / "genome.fa.gz",
    output:
        MARK_DUPLICATES / "{sample}.cram",
    log:
        MARK_DUPLICATES / "{sample}.cram.log",
    threads: 24
    conda:
        "__environment__.yml"
    shell:
        """
        samtools view \
            --threads {threads} \
            --output-fmt CRAM \
            --reference {input.reference} \
            --output {output} \
            {input.bam} \
        2> {log} 1>&2
        """


rule align__mark_duplicates__all:
    """Mark duplicates in all chromosomes and all libraries"""
    input:
        [MARK_DUPLICATES / f"{sample}.cram" for sample in SAMPLES],
