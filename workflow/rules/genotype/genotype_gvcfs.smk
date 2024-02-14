rule genotype__genotype_gvcfs__:
    """Genotype a single chromosome"""
    input:
        vcf_gz=COMBINE_GVCFS / "{chromosome}.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf_gz=GENOTYPE_GVCFS / "{chromosome}.vcf.gz",
    log:
        GENOTYPE_GVCFS / "{chromosome}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["gatk4"]["genotype_gvcfs"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk GenotypeGVCFs \
            {params.extra} \
            --variant {input.vcf_gz} \
            --reference {input.reference} \
            --output {output.vcf_gz} \
        2> {log} 1>&2
        """


rule genotype__genotype_gvcfs__all:
    """Genotype all chromosomes"""
    input:
        [GENOTYPE_GVCFS / f"{chromosome}.vcf.gz" for chromosome in CHROMOSOMES],
