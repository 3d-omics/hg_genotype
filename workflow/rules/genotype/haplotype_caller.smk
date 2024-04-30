rule genotype__haplotype_caller__:
    """Call variants for a single library and chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        cram=RECALIBRATE / "{sample_id}.cram",
        crai=RECALIBRATE / "{sample_id}.cram.crai",
        dict_=REFERENCE / "genome.dict",
    output:
        gvcf_gz=HAPLOTYPE_CALLER / "{sample_id}" / "{region}.gvcf.gz",
    log:
        HAPLOTYPE_CALLER / "{sample_id}" / "{region}.log",
    conda:
        "__environment__.yml"
    params:
        extra=params["gatk4"]["haplotype_caller"]["extra"],
        ploidy=get_ploidy_of_sample_and_chromosome,
        interval=get_interval_for_haplotype_caller,
    group:
        "genotype"
    shell:
        """
        gatk HaplotypeCaller \
            {params.extra} \
            --reference {input.reference} \
            --input {input.cram} \
            --output {output.gvcf_gz} \
            --emit-ref-confidence GVCF \
            --intervals {params.interval} \
            -ploidy {params.ploidy} \
        2> {log} 1>&2
        """


rule genotype__haplotype_caller:
    """Call variants for all libraries and chromosomes"""
    input:
        [
            HAPLOTYPE_CALLER / f"{sample_id}" / f"{region}.gvcf.gz"
            for sample_id in SAMPLES
            for region in REGIONS
        ],
