rule report__sample__:
    """Make a MultiQC report for a single library"""
    input:
        get_reads_fastqc_for_library_report,
        get_reports_from_bwa_for_library_report,
        get_picard_markduplicates_for_library_report,
        RECALIBRATE / "{sample}.bsqr.txt",
        get_gatk4_apply_bqsr_for_library_report,
    output:
        SAMPLE / "{sample}.html",
    log:
        SAMPLE / "{sample}.log",
    conda:
        "__environment__.yml"
    params:
        title=lambda w: f"{w.sample}",
        out_dir=SAMPLE,
    shell:
        """
        multiqc \
            --title {params.title} \
            --force \
            --dirs \
            --dirs-depth 1 \
            --filename {params.title} \
            --outdir {params.out_dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__sample__all:
    """Make a MultiQC report for every library"""
    input:
        [SAMPLE / f"{sample}.html" for sample in SAMPLES],


rule report__sample:
    """Make all MultiQC reports per library"""
    input:
        rules.report__sample__all.input,
