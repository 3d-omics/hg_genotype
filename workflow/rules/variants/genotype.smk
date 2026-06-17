rule variants__genotype__genotype_gvcfs:
    """Genotype a single region"""
    input:
        gvcf=CALL / "{region}.vcf.gz",
        ref=REFERENCE / f"{HOST_NAME}.fa.gz",
    output:
        vcf=GENOTYPE / "{region}.vcf.gz",
    log:
        GENOTYPE / "{region}.log",
    retries: 5
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=7 * 24 * 60,
    wrapper:
        "v7.9.1/bio/gatk/genotypegvcfs"


rule variants__genotype__genotype_gvcfs__all:
    input:
        [GENOTYPE / f"{region}.vcf.gz" for region in REGIONS],


rule variants__genotype__merge_vcfs:
    """Join all the GVCFs into a single one

    Mysteriously MergeVcfs fucks up the file
    """
    input:
        calls=[GENOTYPE / f"{region}.vcf.gz" for region in REGIONS],
    output:
        vcf_gz=GENOTYPE / "all.vcf.gz",
    log:
        GENOTYPE / "all.log",
    params:
        extra="--allow-overlaps",
    wrapper:
        "v7.9.1/bio/bcftools/concat"


rule variants__genotype__merge_vcfs__all:
    input:
        vcf_gz=GENOTYPE / "all.vcf.gz",


rule variants__genotype__all:
    input:
        rules.variants__genotype__genotype_gvcfs__all.input,
        rules.variants__genotype__merge_vcfs__all.input,
