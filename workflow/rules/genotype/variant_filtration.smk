rule genotype__variant_filtration__:
    """Filter variants for a single chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        vcf=get_input_vcf_for_genotype__variant_filtration,
    output:
        vcf=VARIANT_FILTRATION / "{region}.vcf.gz",
    log:
        VARIANT_FILTRATION / "{region}.log",
    conda:
        "__environment__.yml"
    params:
        filter_name=params["gatk4"]["variant_filtration"]["filter_name"],
        filter_expression=params["gatk4"]["variant_filtration"]["filter_expression"],
        extra=params["gatk4"]["variant_filtration"]["extra"],
    group:
        "genotype"
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


rule genotype__variant_filtration__merge__:
    """Merge all VCF chromosomes"""
    input:
        expand(
            VARIANT_FILTRATION / "{region}.vcf.gz",
            region=DIPLOID_REGIONS,
        ),
    output:
        VARIANT_FILTRATION / "variants_filtered.vcf.gz",
    log:
        VARIANT_FILTRATION / "variants_filtered.log",
    conda:
        "__environment__.yml"
    shell:
        """
        bcftools concat \
            --output {output} \
            --output-type z \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """


rule genotype__variant_filtration:
    input:
        rules.genotype__variant_filtration__merge__.output,
