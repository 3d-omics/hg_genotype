rule align__mark_duplicates__:
    """Mark duplicates in a single chromosome from a single library"""
    input:
        bam=MAP / "{sample}.{library}.cram",
        reference=REFERENCE / "genome.fa.gz",
    output:
        bam=MARK_DUPLICATES / "{sample}.{library}.bam",
        metrics=MARK_DUPLICATES / "{sample}.{library}.metrics.tsv",
    log:
        MARK_DUPLICATES / "{sample}.{library}.log",
    conda:
        "__environment__.yml"
    resources:
        mem_mb=8000,
        runtime=360,
    shell:
        """
        gatk MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            --COMPRESSION_LEVEL 1 \
            --REFERENCE_SEQUENCE {input.reference} \
        2> {log} 1>&2
        """


rule align__mark_duplicates__all:
    """Mark duplicates in all chromosomes and all libraries"""
    input:
        [
            MARK_DUPLICATES / f"{sample}.{library}.bam"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
        ],
