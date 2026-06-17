use rule helpers__multiqc as swaps__multiqc with:
    input:
        SOMALIER / "relate.pairs.tsv",
        SOMALIER / "relate.samples.tsv",
    output:
        html=RESULTS / "swaps.html",
        zip=RESULTS / "swaps.zip",
    log:
        RESULTS / "swaps.log",
    params:
        extra="--title swaps --dirs --fullnames --fn_as_s_name --force",


rule swaps__multiqc__all:
    """Collect all per step reports for the pipeline"""
    input:
        RESULTS / "swaps.html",
