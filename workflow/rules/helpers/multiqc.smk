rule helpers__multiqc:
    """Base rule for collecting per-step QC reports with MultiQC."""
    output:
        html=RESULTS / "multiqc.html",
        zip=RESULTS / "multiqc.zip",
    log:
        RESULTS / "multiqc.log",
    resources:
        mem_mb=8 * 1024,
        runtime=2 * 60,
    params:
        extra="--dirs --fullnames --fn_as_s_name --force",
    wrapper:
        "v7.9.1/bio/multiqc"
