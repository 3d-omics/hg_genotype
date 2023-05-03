rule report_one_sample:
    input:
        READS / "{sample}.{library}_1_fastqc.zip",
        READS / "{sample}.{library}_2_fastqc.zip",
        FASTP / "{sample}.{library}_fastp.json",
        BOWTIE2 / "{sample}.{library}.stats.tsv",
        BOWTIE2 / "{sample}.{library}.flagstats.txt",
        BOWTIE2 / "{sample}.{library}.idxstats.tsv",
        PICARD / "{sample}.{library}.stats.tsv",
        PICARD / "{sample}.{library}.flagstats.txt",
        PICARD / "{sample}.{library}.idxstats.tsv",
        PICARD / "{sample}.{library}.metrics.tsv",
        GATK / "{sample}.{library}.base_recalibrator.txt",
    output:
        REPORTS_BY_SAMPLE / "{sample}.{library}.html",
    log:
        REPORTS_BY_SAMPLE / "{sample}.{library}.log",
    conda:
        "../envs/report.yml"
    params:
        library="{sample}.{library}",
        out_dir=REPORTS_BY_SAMPLE,
    shell:
        """
        multiqc \
            --filename {params.library} \
            --outdir {params.out_dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_samples:
    input:
        [
            REPORTS_BY_SAMPLE / f"{sample}.{library}.html"
            for sample, library in SAMPLE_LIB
        ],
