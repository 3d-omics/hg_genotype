rule helpers__bcftools__index_vcf_gz__:
    """Index a vcf.gz"""
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.csi",
    conda:
        "__environment__.yml"
    log:
        "{prefix}.vcf.gz.csi.log",
    shell:
        "bcftools index {input} 2> {log} 1>&2"
