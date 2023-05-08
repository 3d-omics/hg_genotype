def get_picard_per_sample_files(wildcards):
    sample = wildcards.sample
    library = wildcards.library
    ANALYSES = ["stats.tsv", "flagstats.txt", "idxstats.tsv", "metrics.tsv"]
    files = [
        PICARD / f"markduplicates/{sample}.{library}.{chromosome}.{analysis}"
        for chromosome in CHROMOSOMES
        for analysis in ANALYSES
    ]
    return files


def get_gatk4_base_recalibrator_per_sample_files(wildcards):
    files = [
        GATK
        / f"base_recalibrator/{wildcards.sample}.{wildcards.library}.{chromosome}.txt"
        for chromosome in CHROMOSOMES
    ]
    return files


rule report_sample_one:
    input:
        READS / "{sample}.{library}_1_fastqc.zip",
        READS / "{sample}.{library}_2_fastqc.zip",
        FASTP / "{sample}.{library}_fastp.json",
        BOWTIE2 / "{sample}.{library}.stats.tsv",
        BOWTIE2 / "{sample}.{library}.flagstats.txt",
        BOWTIE2 / "{sample}.{library}.idxstats.tsv",
        get_picard_per_sample_files,
        get_gatk4_base_recalibrator_per_sample_files,
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
            --title {params.library} \
            --force \
            --filename {params.library} \
            --outdir {params.out_dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_sample_all:
    input:
        [
            REPORTS_BY_SAMPLE / f"{sample}.{library}.html"
            for sample, library in SAMPLE_LIB
        ],


rule report_samples:
    input:
        rules.report_sample_all.input,
