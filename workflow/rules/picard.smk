rule picard_markduplicates:
    input:
        cram=BOWTIE2 / "{sample}.{library}.cram",
        reference=REFERENCE / "genome.fa.gz",
    output:
        bam=PICARD / "{sample}.{library}.bam",
        metrics=PICARD / "{sample}.{library}.metrics.tsv",
    log:
        PICARD / "{sample}.{library}.log",
    conda:
        "../envs/picard.yml"
    threads: 1  # TODO: tune it
    params:
        None,
    shell:
        """
        picard MarkDuplicates \
            --INPUT {input.cram} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            --COMPRESSION_LEVEL 1 \
            --REFERENCE_SEQUENCE {input.reference} \
        2> {log} 1>&2
        """


rule picard:
    input:
        [PICARD / f"{sample}.{library}.bam" for sample, library in SAMPLE_LIB],


rule picard_reports:
    input:
        [
            PICARD / f"{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in "stats.tsv flagstats.txt idxstats.tsv metrics.tsv".split()
        ],
