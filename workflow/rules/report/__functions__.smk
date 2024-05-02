def get_reads_fastqc_for_library_report(wildcards):
    sample_id = wildcards.sample_id
    libraries = samples[samples["sample_id"] == sample_id]["library_id"].values
    files = [
        READS / f"{sample_id}.{library_id}_{end}_fastqc.html"
        for library_id in libraries
        for end in ["1", "2"]
    ]
    return files


def get_reports_from_bwa_for_library_report(wildcards):
    """Get all bwa reports for a single library"""
    sample_id = wildcards.sample_id
    libraries = samples[samples["sample_id"] == sample_id]["library_id"].values
    files = [
        MAP / f"{sample_id}.{library_id}.{report}"
        for report in BAM_REPORTS
        for library_id in libraries
    ]
    return files


def get_picard_markduplicates_for_library_report(wildcards):
    """Get all picard markduplicates reports for a single library"""
    files = [
        MARK_DUPLICATES / f"{wildcards.sample_id}.{report}" for report in BAM_REPORTS
    ]
    return files


def get_gatk4_apply_bqsr_for_library_report(wildcards):
    """Get gatk4 apply bqsr reports for a single library"""
    return [RECALIBRATE / f"{wildcards.sample_id}.{report}" for report in BAM_REPORTS]
