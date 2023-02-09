rule bcftools_index_vcf_gz:
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.csi",
    conda:
        "../envs/bcftools.yml"
    log:
        "{prefix}.vcf.gz.csi.log",
    shell:
        "bcftools index {input} 2> {log} 1>&2"
