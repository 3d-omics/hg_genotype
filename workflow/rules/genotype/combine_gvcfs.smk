rule genotype__combine_gvcfs__:
    """Combine gVCFs to get a region"""
    input:
        vcf_gzs=get_files_to_genotype,
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf_gz=COMBINE_GVCFS / "{region}.vcf.gz",
    log:
        COMBINE_GVCFS / "{region}.log",
    conda:
        "__environment__.yml"
    params:
        variant_line=compose_v_line,
        extra=params["gatk4"]["combine_gvcfs"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    group:
        "genotype"
    shell:
        """
        gatk CombineGVCFs \
            {params.extra} \
            {params.variant_line} \
            --reference {input.reference} \
            --output {output.vcf_gz} \
        2> {log} 1>&2
        """


rule genotype__combine_gvcfs:
    """Get all chromosomal gVCFs"""
    input:
        [COMBINE_GVCFS / f"{region}.vcf.gz" for region in REGIONS],
