rule report__step__reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads__fastqc__all.input,
    output:
        html=STEP / "reads.html",
    log:
        STEP / "reads.log",
    conda:
        "__environment__.yml"
    params:
        dir=STEP,
    shell:
        """
        multiqc \
            --filename reads \
            --title reads \
            --force \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__align:
    """Collect all reports for the bowtie2 step"""
    input:
        [
            MAP / f"{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
        ],
        [
            MARK_DUPLICATES / f"{sample}.{library}" / f"{chromosome}.{report}"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
            for report in PICARD_REPORTS
        ],
    output:
        html=STEP / "align.html",
    log:
        STEP / "align.log",
    conda:
        "__environment__.yml"
    params:
        dir=STEP,
    shell:
        """
        multiqc \
            --title align \
            --force \
            --filename align \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__genotype:
    """Collect all reports for the gatk4 step"""
    input:
        rules.gatk4_report.input,
    output:
        html=STEP / "genotype.html",
    log:
        STEP / "genotype.log",
    conda:
        "__environment__.yml"
    params:
        dir=STEP,
    shell:
        """
        multiqc \
            --title genotype \
            --force \
            --filename genotype \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report__step__annotate:
    """Collect all reports for the snpeff step"""
    input:
        rules.annotate__snpeff.input,
    output:
        html=STEP / "annotate.html",
    log:
        STEP / "annotate.log",
    conda:
        "__environment__.yml"
    params:
        dir=STEP,
    shell:
        """
        multiqc \
            --title annotate \
            --force \
            --filename annotate \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report__step__reads.output,
        rules.report__step__align.output,
        rules.report__step__genotype.output,
        rules.report__step__annotate.output,
