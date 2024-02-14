rule genotype__combine_gvcfs__:
    """Combine gVCFs to get a chromosome"""
    input:
        vcf_gzs=get_files_to_genotype,
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf_gz=COMBINE_GVCFS / "{chromosome}.vcf.gz",
    log:
        COMBINE_GVCFS / "{chromosome}.log",
    conda:
        "__environment__.yml"
    params:
        variant_line=compose_v_line,
        extra=params["gatk4"]["combine_gvcfs"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk CombineGVCFs \
            {params.extra} \
            {params.variant_line} \
            --reference {input.reference} \
            --output {output.vcf_gz} \
        2> {log} 1>&2
        """


rule genotype__combine_gvcfs__all:
    """Get all chromosomal gVCFs"""
    input:
        [COMBINE_GVCFS / f"{chromosome}.vcf.gz" for chromosome in CHROMOSOMES],
