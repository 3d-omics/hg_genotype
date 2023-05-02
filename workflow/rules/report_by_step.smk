rule reports_step_reads:
    input:
        [
            READS / f"{sample}.{library}_{end}_fastqc.zip"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
        ],
    output:
        html=REPORTS_BY_STEP / "raw.html",
    log:
        REPORTS_BY_STEP / "raw.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORTS_BY_STEP,
    shell:
        """
        multiqc \
            --filename raw \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule reports_step_fastp:
    input:
        rules.fastp_fastqc.input,
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
            --filename fastp \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule reports_step_bowtie2:
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
            --filename bowtie2 \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule reports_step:
    input:
        rules.reports_step_reads.output,
        rules.reports_step_fastp.output,
        rules.reports_step_bowtie2.output,