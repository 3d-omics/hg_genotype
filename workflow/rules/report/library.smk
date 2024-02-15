rule report__library__:
    """Make a MultiQC report for a single library"""
    input:
        READS / "{sample}.{library}_1_fastqc.zip",
        READS / "{sample}.{library}_2_fastqc.zip",
        MAP / "{sample}.{library}.stats.tsv",
        MAP / "{sample}.{library}.flagstats.txt",
        MAP / "{sample}.{library}.idxstats.tsv",
        get_picard_markduplicates_for_library_report,
        get_gatk4_base_recalibrator_for_library_report,
        get_gatk4_apply_bqsr_for_library_report,
    output:
        LIBRARY / "{sample}.{library}.html",
    log:
        LIBRARY / "{sample}.{library}.log",
    conda:
        "__environment__.yml"
    params:
        library="{sample}.{library}",
        out_dir=LIBRARY,
    shell:
        """
        multiqc \
            --title {params.library} \
            --force \
            --filename {params.library} \
            --outdir {params.out_dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__library__all:
    """Make a MultiQC report for every library"""
    input:
        [LIBRARY / f"{sample}.{library}.html" for sample, library in SAMPLE_LIB],


rule report__library:
    """Make all MultiQC reports per library"""
    input:
        rules.report__library__all.input,
