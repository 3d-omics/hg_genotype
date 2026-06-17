use rule helpers__multiqc as align__multiqc with:
    input:
        reads=[
            READS / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
        bwamem2=[
            MAP / f"{sample_id}.{library_id}.stats"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        markduplicates=[MARK_DUPLICATES / f"{sample_id}.stats" for sample_id in SAMPLES],
        recalibrate=[
            RECALIBRATE / f"{sample_id}.stats"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        html=RESULTS / "align.html",
        zip=RESULTS / "align.zip",
    log:
        RESULTS / "align.log",
    params:
        extra="--title align --dirs --fullnames --fn_as_s_name --force",


rule align__multiqc__all:
    """Collect all reports for the align submodule step"""
    input:
        RESULTS / "align.html",
