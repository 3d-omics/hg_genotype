rule align__recalibrate__baserecalibrator:
    """Compute the recalibration table for a single library and chromosome"""
    input:
        cram=MARK_DUPLICATES / "{sample}.{library}.cram",
        crai=MARK_DUPLICATES / "{sample}.{library}.cram.crai",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        known_sites=REFERENCE / "known_variants.vcf.gz",
        tbi=REFERENCE / "known_variants.vcf.gz.tbi",
    output:
        table=RECALIBRATE / "{sample}.{library}.bsqr.txt",
    log:
        RECALIBRATE / "{sample}.{library}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["gatk4"]["base_recalibrator"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk BaseRecalibrator \
            {params.extra} \
            --input {input.cram} \
            --reference {input.reference} \
            --known-sites {input.known_sites} \
            --output {output.table} \
        2> {log} 1>&2
        """


rule align__recalibrate__applybqsr:
    """Apply the recalibration table to a single library and chromosome"""
    input:
        cram=MARK_DUPLICATES / "{sample}.{library}.cram",
        reference=REFERENCE / "genome.fa.gz",
        table=RECALIBRATE / "{sample}.{library}.bsqr.txt",
        dict_=REFERENCE / "genome.dict",
    output:
        bam=pipe(RECALIBRATE / "{sample}.{library}.bam"),
    log:
        RECALIBRATE / "{sample}.{library}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["gatk4"]["apply_bqsr"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    threads: 0
    shell:
        """
        gatk ApplyBQSR \
            {params.extra} \
            --input {input.cram} \
            --reference {input.reference} \
            --bqsr-recal-file {input.table} \
            --output {output.bam} \
        2> {log} 1>&2
        """


rule align__recalibrate__bam_to_cram:
    input:
        bam=RECALIBRATE / "{sample}.{library}.bam",
        reference=REFERENCE / "genome.fa.gz",
    output:
        cram=RECALIBRATE / "{sample}.{library}.cram",
    log:
        RECALIBRATE / "{sample}.{library}.cram.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        memory="4G",
        runtime=1440,
    shell:
        """
        samtools view \
            --threads {threads} \
            --output-fmt CRAM \
            --reference {input.reference} \
            --output {output.cram} \
            {input.bam} \
        2> {log} 1>&2
        """


rule align__recalibrate__all:
    """Compute recalibration for all chromosomes and libraries"""
    input:
        [
            RECALIBRATE / f"{sample}.{library}.cram"
            for sample, library in SAMPLE_LIB
            for chromosome in CHROMOSOMES
        ],
