def get_picard_markduplicates_for_chromosome_report(wildcards):
    """Get all picard markduplicates reports for a single chromosome"""
    chromosome = wildcards.chromosome
    files = [
        PICARD / f"markduplicates/{sample}.{library}.{chromosome}.{report}"
        for sample, library in SAMPLE_LIB
        for report in PICARD_REPORTS
    ]
    return files


def get_gatk4_base_recalibrator_for_chromosome_report(wildcards):
    """Get gatk4 base recalibrator reports for a single chromosome"""
    chromosome = wildcards.chromosome
    files = [
        GATK / f"base_recalibrator/{sample}.{library}.{chromosome}.txt"
        for sample, library in SAMPLE_LIB
    ]
    return files


def get_gatk4_apply_bqsr_for_chromosome_report(wildcards):
    """Get gatk4 apply bqsr reports for a single chromosome"""
    chromosome = wildcards.chromosome
    files = [
        GATK / f"apply_bqsr/{sample}.{library}.{chromosome}.{report}"
        for sample, library in SAMPLE_LIB
        for report in BAM_REPORTS
    ]
    return files
