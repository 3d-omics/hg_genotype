rule genotype__calculate_genotype_posteriors__:
    """Calculate genotype posteriors for a single chromosome"""
    input:
        vcf=GENOTYPE_GVCFS / "{region}.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf=POSTERIORS / "{region}.vcf.gz",
    log:
        POSTERIORS / "{region}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["gatk4"]["calculate_genotype_posteriors"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    group:
        "genotype"
    shell:
        """
        gatk CalculateGenotypePosteriors \
            {params.extra} \
            --output {output.vcf} \
            --variant {input.vcf} \
            --reference {input.reference} \
        2> {log} 1>&2
        """


rule genotype__calculate_genotype_posteriors__all:
    """Calculate genotype posteriors for all chromosomes"""
    input:
        [POSTERIORS / f"{region}.vcf.gz" for region in DIPLOID_REGIONS],
