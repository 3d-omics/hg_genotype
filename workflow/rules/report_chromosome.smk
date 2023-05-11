rule report_chromosome_one:
    """Generate a report for a single chromosome"""
    input:
        get_picard_markduplicates_per_chromosome_files,
        get_gatk4_base_recalibrator_per_chromosome_files,
    output:
        REPORT_CHR / "{chromosome}.html",
    log:
        REPORT_CHR / "{chromosome}.log",
    conda:
        "../envs/report.yml"
    params:
        chromosome="{chromosome}",
        out_dir=REPORT_CHR,
    shell:
        """
        multiqc \
            --title {params.chromosome} \
            --force \
            --filename {params.chromosome} \
            --outdir {params.out_dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_chromosome_all:
    """Generate all reports for all chromosomes"""
    input:
        [REPORT_CHR / f"{chromosome}.html" for chromosome in CHROMOSOMES],


rule report_chromosome:
    """Generate all report for all chromosomes"""
    input:
        rules.report_chromosome_all.input,
