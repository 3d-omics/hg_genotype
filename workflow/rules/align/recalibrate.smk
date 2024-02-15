rule align__recalibrate__baserecalibrator:
    """Compute the recalibration table for a single library and chromosome"""
    input:
        bam=MARK_DUPLICATES / "{sample}.{library}.bam",
        bai=MARK_DUPLICATES / "{sample}.{library}.bam.bai",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        known_sites=REFERENCE / "known_variants.vcf.gz",
        csi=REFERENCE / "known_variants.vcf.gz.tbi",
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
            --input {input.bam} \
            --reference {input.reference} \
            --known-sites {input.known_sites} \
            --output {output.table} \
        2> {log} 1>&2
        """


rule align__recalibrate__baserecalibrator__all:
    """Compute recalibration for all chromosomes and libraries"""
    input:
        [
            RECALIBRATE / f"{sample}.{library}.bsqr.txt"
            for sample, library in SAMPLE_LIB
            for chromosome in CHROMOSOMES
        ],


rule align__recalibrate__applybqsr:
    """Apply the recalibration table to a single library and chromosome"""
    input:
        bam=MARK_DUPLICATES / "{sample}.{library}.bam",
        reference=REFERENCE / "genome.fa.gz",
        table=RECALIBRATE / "{sample}.{library}.bsqr.txt",
        dict_=REFERENCE / "genome.dict",
    output:
        bam=RECALIBRATE / "{sample}.{library}.bam",
    log:
        RECALIBRATE / "{sample}.{library}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["gatk4"]["apply_bqsr"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk ApplyBQSR \
            {params.extra} \
            --input {input.bam} \
            --reference {input.reference} \
            --bqsr-recal-file {input.table} \
            --output {output.bam} \
        2> {log} 1>&2
        """


rule align__recalibrate__applybqsr__all:
    """Compute recalibration for all chromosomes and libraries"""
    input:
        [
            RECALIBRATE / f"{sample}.{library}.bam"
            for sample, library in SAMPLE_LIB
            for chromosome in CHROMOSOMES
        ],
