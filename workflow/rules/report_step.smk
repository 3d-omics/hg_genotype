rule report_step_reads:
    input:
        rules.reads_fastqc_all.input,
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "../envs/report.yml"
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


rule report_step_fastp:
    input:
        rules.fastp_report_all.input,
    output:
        html=REPORT_STEP / "fastp.html",
    log:
        REPORT_STEP / "fastp.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
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
        rules.bowtie2_report_all.input,
    output:
        html=REPORT_STEP / "bowtie2.html",
    log:
        REPORT_STEP / "bowtie2.log",
    conda:
        "../envs/report.yml"
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


rule report_step_somalier:
    input:
        rules.somalier_report.input,
    output:
        html=REPORT_STEP / "somalier.html",
    log:
        REPORT_STEP / "somalier.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title somalier \
            --force \
            --filename somalier \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_picard:
    input:
        rules.picard_report_all.input,
    output:
        html=REPORT_STEP / "picard.html",
    log:
        REPORT_STEP / "picard.log",
    conda:
        "../envs/report.yml"
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
    input:
        rules.gatk4_base_recalibrator_all.input,
    output:
        html=REPORT_STEP / "gatk4.html",
    log:
        REPORT_STEP / "gatk4.log",
    conda:
        "../envs/report.yml"
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
    input:
        rules.snpeff_report.input,
    output:
        html=REPORT_STEP / "snpeff.html",
    log:
        REPORT_STEP / "snpeff.log",
    conda:
        "../envs/report.yml"
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
    input:
        rules.report_step_reads.output,
        rules.report_step_fastp.output,
        rules.report_step_bowtie2.output,
        rules.report_step_somalier.output,
        rules.report_step_picard.output,
        rules.report_step_gatk4.output,
        rules.report_step_snpeff.output,
