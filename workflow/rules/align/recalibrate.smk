rule align__recalibrate__baserecalibrator:
    """Compute the recalibration table for a single sample"""
    input:
        bam=MARK_DUPLICATES / "{sample_id}.cram",
        crai=MARK_DUPLICATES / "{sample_id}.cram.crai",
        ref=ancient(REFERENCE / f"{HOST_NAME}.fa.gz"),
        dict=REFERENCE / f"{HOST_NAME}.dict",
        known=ancient(REFERENCE / f"{HOST_NAME}.vcf.gz"),
        tbi=REFERENCE / f"{HOST_NAME}.vcf.gz.tbi",
    output:
        recal_table=temp(RECALIBRATE / "{sample_id}.bsqr.txt"),
    log:
        RECALIBRATE / "{sample_id}.baserecalibrator.log",
    benchmark:
        RECALIBRATE / "{sample_id}.baserecalibrator.benchmark.tsv"
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
    wrapper:
        "v7.9.1/bio/gatk/baserecalibrator"


rule align__recalibrate__applybqsr:
    """Apply the recalibration table to a single sample"""
    input:
        bam=MARK_DUPLICATES / "{sample_id}.cram",
        ref=ancient(REFERENCE / f"{HOST_NAME}.fa.gz"),
        dict=REFERENCE / f"{HOST_NAME}.dict",
        recal_table=RECALIBRATE / "{sample_id}.bsqr.txt",
    output:
        bam=RECALIBRATE / "{sample_id}.cram",
    log:
        RECALIBRATE / "{sample_id}.applybqsr.log",
    benchmark:
        RECALIBRATE / "{sample_id}.applybqsr.benchmark.tsv"
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
    wrapper:
        "v7.9.1/bio/gatk/applybqsr"


rule align__recalibrate__all:
    """Compute recalibration for all samples"""
    input:
        [RECALIBRATE / f"{sample_id}.cram" for sample_id in SAMPLES],
