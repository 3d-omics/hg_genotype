rule genotype__variant_filtration__:
    """Filter variants for a single chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        vcf=GATK / "variants_posteriors/{chromosome}.vcf.gz",
    output:
        vcf=GATK / "variants_filtered/{chromosome}.vcf.gz",
    log:
        GATK / "variants_filtered/{chromosome}.log",
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
            GATK / "variants_filtered/{chromosome}.vcf.gz",
            chromosome=CHROMOSOMES,
        ),
    output:
        GATK / "variants_filtered.vcf.gz",
    log:
        GATK / "variants_filtered.log",
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
