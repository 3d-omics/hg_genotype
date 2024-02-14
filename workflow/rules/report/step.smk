rule report_step_reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads__fastqc__all.input,
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
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


rule report_step_bowtie2:
    """Collect all reports for the bowtie2 step"""
    input:
        [
            MAP / f"{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
        ],
    output:
        html=REPORT_STEP / "bowtie2.html",
    log:
        REPORT_STEP / "bowtie2.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title bowtie2 \
            --force \
            --filename bowtie2 \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_picard:
    """Collect all reports for the picard step"""
    input:
        [
            MARK_DUPLICATES / f"{sample}.{library}" / f"{chromosome}.{report}"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
            for report in PICARD_REPORTS
        ],
    output:
        html=REPORT_STEP / "picard.html",
    log:
        REPORT_STEP / "picard.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title picard \
            --force \
            --filename picard \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_gatk4:
    """Collect all reports for the gatk4 step"""
    input:
        rules.gatk4_report.input,
    output:
        html=REPORT_STEP / "gatk4.html",
    log:
        REPORT_STEP / "gatk4.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title gatk4 \
            --force \
            --filename gatk4 \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_snpeff:
    """Collect all reports for the snpeff step"""
    input:
        rules.snpeff_report.input,
    output:
        html=REPORT_STEP / "snpeff.html",
    log:
        REPORT_STEP / "snpeff.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title snpeff \
            --force \
            --filename snpeff \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report_step_reads.output,
        # rules.report_step_fastp.output,
        rules.report_step_bowtie2.output,
        rules.report_step_picard.output,
        rules.report_step_gatk4.output,
        rules.report_step_snpeff.output,
