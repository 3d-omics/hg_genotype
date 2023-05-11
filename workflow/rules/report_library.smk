rule report_library_one:
    input:
        READS / "{sample}.{library}_1_fastqc.zip",
        READS / "{sample}.{library}_2_fastqc.zip",
        FASTP / "{sample}.{library}_fastp.json",
        BOWTIE2 / "{sample}.{library}.stats.tsv",
        BOWTIE2 / "{sample}.{library}.flagstats.txt",
        BOWTIE2 / "{sample}.{library}.idxstats.tsv",
    output:
        REPORT_LIBRARY / "{sample}.{library}.html",
    log:
        REPORT_LIBRARY / "{sample}.{library}.log",
    conda:
        "../envs/report.yml"
    params:
        library="{sample}.{library}",
        out_dir=REPORT_LIBRARY,
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


rule report_library_all:
    input:
        [REPORT_LIBRARY / f"{sample}.{library}.html" for sample, library in SAMPLE_LIB],


rule report_library:
    input:
        rules.report_library_all.input,
