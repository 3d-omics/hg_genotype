rule report_step_reads:
    input:
        rules.reads_fastqc.input,
    output:
        html=REPORTS_BY_STEP / "reads.html",
    log:
        REPORTS_BY_STEP / "reads.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORTS_BY_STEP,
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


rule report_step_fastp:
    input:
        rules.fastp_reports.input,
    output:
        html=REPORTS_BY_STEP / "fastp.html",
    log:
        REPORTS_BY_STEP / "fastp.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORTS_BY_STEP,
    shell:
        """
        multiqc \
            --title fastp \
            --force \
            --filename fastp \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_bowtie2:
    input:
        rules.bowtie2_reports.input,
    output:
        html=REPORTS_BY_STEP / "bowtie2.html",
    log:
        REPORTS_BY_STEP / "bowtie2.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORTS_BY_STEP,
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
    input:
        rules.picard_reports.input,
    output:
        html=REPORTS_BY_STEP / "picard.html",
    log:
        REPORTS_BY_STEP / "picard.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORTS_BY_STEP,
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
    input:
        rules.gatk4_base_recalibrator_all.input,
    output:
        html=REPORTS_BY_STEP / "gatk4.html",
    log:
        REPORTS_BY_STEP / "gatk4.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORTS_BY_STEP,
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
    input:
        rules.snpeff_report.input,
    output:
        html=REPORTS_BY_STEP / "snpeff.html",
    log:
        REPORTS_BY_STEP / "snpeff.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORTS_BY_STEP,
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


rule report_steps:
    input:
        rules.report_step_reads.output,
        rules.report_step_fastp.output,
        rules.report_step_bowtie2.output,
        rules.report_step_picard.output,
        rules.report_step_gatk4.output,
        rules.report_step_snpeff.output,
