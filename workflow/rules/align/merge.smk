rule align__merge__:
    input:
        bams=lambda w: [
            RECALIBRATE / f"{w.sample}.{w.library}" / f"{chromosome}.bam"
            for chromosome in CHROMOSOMES
        ],
        reference=REFERENCE / "genome.fa.gz",
    output:
        cram=MERGE / "{sample}.{library}.cram",
    log:
        MERGE / "{sample}.{library}.log",
    conda:
        "__environment__.yml"
    resources:
        mem_mb=8000,
        runtime=1440,
    threads: 8
    shell:
        """
        samtools merge \
            --output-fmt CRAM \
            --reference {input.reference} \
            --threads {threads} \
            -l 9 \
            -o {output.cram} \
            {input.bams} \
        2> {log} 1>&2
        """


rule align__merge__all:
    """Apply the recalibration table to all libraries and chromosomes"""
    input:
        [
            MERGE / f"{sample}.{library}.cram"
            for sample, library in SAMPLE_LIB
            for chromosome in CHROMOSOMES
        ],
