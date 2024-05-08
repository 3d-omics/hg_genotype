rule report__sample__:
    """Make a MultiQC report for a single library"""
    input:
        get_reads_fastqc_for_library_report,
        get_reports_from_bwa_for_library_report,
        get_picard_markduplicates_for_library_report,
        RECALIBRATE / "{sample_id}.bsqr.txt",
        get_gatk4_apply_bqsr_for_library_report,
        VEP / "{sample_id}.vep.html",
    output:
        SAMPLE / "{sample_id}.html",
    log:
        SAMPLE / "{sample_id}.log",
    conda:
        "__environment__.yml"
    params:
        title=lambda w: f"{w.sample_id}",
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
        [SAMPLE / f"{sample_id}.html" for sample_id in SAMPLES],


rule report__sample:
    """Make all MultiQC reports per library"""
    input:
        rules.report__sample__all.input,
