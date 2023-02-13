rule fastqc:
    input:
        "{prefix}.fq.gz",
    output:
        html="{prefix}_fastqc.html",
        zip="{prefix}_fastqc.zip",
    params:
        "--quiet",
    log:
        "{prefix}_fastqc.log",
    threads: 1
    wrapper:
        "v1.23.1/bio/fastqc"
