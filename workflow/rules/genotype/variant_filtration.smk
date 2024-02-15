def get_input_vcf_for_genotype__variant_filtration(wildcards):
    chromosome = wildcards.chromosome
    return (
        POSTERIORS / f"{chromosome}.vcf.gz"
        if chromosome in DIPLOID_CHROMOSOMES
        else GENOTYPE_GVCFS / f"{chromosome}.vcf.gz"
    )


rule genotype__variant_filtration__:
    """Filter variants for a single chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        vcf=get_input_vcf_for_genotype__variant_filtration,
    output:
        vcf=VARIANT_FILTRATION / "{chromosome}.vcf.gz",
    log:
        VARIANT_FILTRATION / "{chromosome}.log",
    conda:
        "__environment__.yml"
    params:
        filter_name=params["gatk4"]["variant_filtration"]["filter_name"],
        filter_expression=params["gatk4"]["variant_filtration"]["filter_expression"],
        extra=params["gatk4"]["variant_filtration"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk VariantFiltration \
            {params.extra} \
            --reference {input.reference} \
            --variant {input.vcf} \
            --output {output.vcf} \
            --filter-expression '{params.filter_expression}' \
            --filter-name '{params.filter_name}' \
        2> {log} 1>&2
        """


rule genotype__variant_filtration__merge:
    """Merge all VCF chromosomes"""
    input:
        expand(
            VARIANT_FILTRATION / "{chromosome}.vcf.gz",
            chromosome=CHROMOSOMES,
        ),
    output:
        VARIANT_FILTRATION / "variants_filtered.vcf.gz",
    log:
        VARIANT_FILTRATION / "variants_filtered.log",
    conda:
        "__environment__.yml"
    threads: 24
    shell:
        """
        bcftools concat \
            --output {output} \
            --output-type z9 \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """


rule genotype__variant_filtration:
    input:
        rules.genotype__variant_filtration__merge.output,
