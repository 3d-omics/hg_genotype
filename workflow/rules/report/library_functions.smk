def get_picard_markduplicates_for_library_report(wildcards):
    """Get all picard markduplicates reports for a single library"""
    sample = wildcards.sample
    library = wildcards.library
    files = [
        MARK_DUPLICATES / f"{sample}.{library}" / f"{chromosome}.{report}"
        for chromosome in CHROMOSOMES
        for report in PICARD_REPORTS
    ]
    return files


def get_gatk4_base_recalibrator_for_library_report(wildcards):
    """Get gatk4 base recalibrator reports for a single library"""
    sample = wildcards.sample
    library = wildcards.library
    files = [
        RECALIBRATE / f"{sample}.{library}" / f"{chromosome}.bsqr.txt"
        for chromosome in CHROMOSOMES
    ]
    return files


def get_gatk4_apply_bqsr_for_library_report(wildcards):
    """Get gatk4 apply bqsr reports for a single library"""
    sample = wildcards.sample
    library = wildcards.library
    files = [
        RECALIBRATE / f"{sample}.{library}" / f"{chromosome}.{report}"
        for chromosome in CHROMOSOMES
        for report in BAM_REPORTS
    ]
    return files
