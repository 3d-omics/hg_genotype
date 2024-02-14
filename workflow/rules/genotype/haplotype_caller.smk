rule genotype__haplotype_caller__:
    """Call variants for a single library and chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        bam=RECALIBRATE / "{sample}.{library}" / "{chromosome}.bam",
        dict_=REFERENCE / "genome.dict",
    output:
        gvcf_gz=GATK / "haplotype_caller/{sample}.{library}" / "{chromosome}.gvcf.gz",
    log:
        GATK / "haplotype_caller/{sample}.{library}" / "{chromosome}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["gatk4"]["haplotype_caller"]["extra"],
        ploidy=get_ploidy_of_sample_and_chromosome,
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk HaplotypeCaller \
            {params.extra} \
            --reference {input.reference} \
            --input {input.bam} \
            --output {output.gvcf_gz} \
            --emit-ref-confidence GVCF \
            -ploidy {params.ploidy} \
        2> {log} 1>&2
        """


rule genotype__haplotype_caller__all:
    """Call variants for all libraries and chromosomes"""
    input:
        [
            GATK / f"haplotype_caller/{sample}.{library}" / f"{chromosome}.gvcf.gz"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
        ],
