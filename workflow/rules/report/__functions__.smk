def get_reads_fastqc_for_library_report(wildcards):
    sample = wildcards.sample
    libraries = samples[samples["sample"] == sample]["library"].values
    files = [
        READS / f"{sample}.{library}_{end}_fastqc.html"
        for library in libraries
        for end in ["1", "2"]
    ]
    return files


def get_reports_from_bwa_for_library_report(wildcards):
    """Get all bwa reports for a single library"""
    sample = wildcards.sample
    libraries = samples[samples["sample"] == sample]["library"].values
    files = [
        MAP / f"{sample}.{library}.{report}"
        for report in BAM_REPORTS
        for library in libraries
    ]
    return files


def get_picard_markduplicates_for_library_report(wildcards):
    """Get all picard markduplicates reports for a single library"""
    files = [MARK_DUPLICATES / f"{wildcards.sample}.{report}" for report in BAM_REPORTS]
    return files


def get_gatk4_apply_bqsr_for_library_report(wildcards):
    """Get gatk4 apply bqsr reports for a single library"""
    return [RECALIBRATE / f"{wildcards.sample}.{report}" for report in BAM_REPORTS]
