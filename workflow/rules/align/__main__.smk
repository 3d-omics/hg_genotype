include: "__functions__.smk"
include: "index.smk"
include: "map.smk"
include: "split.smk"
include: "mark_duplicates.smk"
include: "recalibrate.smk"


rule picard_report_all:
    """Generate reports for all chromosomes and all libraries"""
    input:
        [
            MARK_DUPLICATES / f"{sample}.{library}" / f"{chromosome}.{report}"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
            for report in PICARD_REPORTS
        ],


rule picard:
    """Run all picard steps and get all reports"""
    input:
        rules.picard_markduplicates_all.input,
        rules.picard_report_all.input,
