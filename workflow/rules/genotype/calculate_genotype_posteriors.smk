rule genotype__calculate_genotype_posteriors__:
    """Calculate genotype posteriors for a single chromosome"""
    input:
        vcf=GATK / "genotyped_variants/{chromosome}.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf=GATK / "variants_posteriors/{chromosome}.vcf.gz",
    log:
        GATK / "variants_posteriors/{chromosome}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["gatk4"]["calculate_genotype_posteriors"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
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
        [
            GATK / f"variants_posteriors/{chromosome}.vcf.gz"
            for chromosome in CHROMOSOMES
        ],
