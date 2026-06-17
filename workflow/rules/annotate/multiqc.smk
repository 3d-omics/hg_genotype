use rule helpers__multiqc as annotate__multiqc with:
    input:
        [VEP / f"{sample}.vep.html" for sample in SAMPLES],
    output:
        html=RESULTS / "annotate.html",
        zip=RESULTS / "annotate.zip",
    log:
        RESULTS / "annotate.log",
    params:
        extra="--title annotate --dirs --fullnames --fn_as_s_name --force",


rule annotate__multiqc__all:
    input:
        RESULTS / "annotate.html",
