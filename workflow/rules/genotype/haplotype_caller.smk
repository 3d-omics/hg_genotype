rule genotype__haplotype_caller__:
    """Call variants for a single library and chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        cram=RECALIBRATE / "{sample}.{library}.cram",
        dict_=REFERENCE / "genome.dict",
    output:
        gvcf_gz=HAPLOTYPE_CALLER / "{sample}.{library}" / "{chromosome}.gvcf.gz",
    log:
        HAPLOTYPE_CALLER / "{sample}.{library}" / "{chromosome}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["gatk4"]["haplotype_caller"]["extra"],
        ploidy=get_ploidy_of_sample_and_chromosome,
        chromosome=lambda w: w.chromosome,
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk HaplotypeCaller \
            {params.extra} \
            --reference {input.reference} \
            --input {input.cram} \
            --output {output.gvcf_gz} \
            --emit-ref-confidence GVCF \
            --intervals {params.chromosome} \
            -ploidy {params.ploidy} \
        2> {log} 1>&2
        """


rule genotype__haplotype_caller__all:
    """Call variants for all libraries and chromosomes"""
    input:
        [
            HAPLOTYPE_CALLER / f"{sample}.{library}" / f"{chromosome}.gvcf.gz"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
        ],
