def get_files_to_genotype(wildcards):
    """Get files to genotype for a sample, library and chromosome"""
    return [
        HAPLOTYPE_CALLER / sample_id / f"{wildcards.region}.gvcf.gz"
        for sample_id in SAMPLES
    ]


def compose_v_line(wildcards):
    """Compose the -v line for gatk4 genotype gvcfs"""
    files = get_files_to_genotype(wildcards)
    text = ""
    for file in files:
        text += f"--variant {file} "
    return text


def get_ploidy_of_sample_and_chromosome(wildcards):
    """Get the ploidy of a sample and chromosome"""
    sample_id = wildcards.sample_id
    region = wildcards.region
    chromosome = REGIONS_BED4[REGIONS_BED4.name == region].chrom.values[0]
    sex = samples[samples["sample_id"] == sample_id][["sex"]]

    male_chrs = features["reference"]["male_chromosomes"]
    female_chrs = features["reference"]["female_chromosomes"]
    mitochondria = features["reference"]["mitochondria"]

    if chromosome in mitochondria:
        return 1
    if (len(male_chrs) == 2) & (chromosome in male_chrs):
        # ZZ won't trigger XY will
        return 1
    if (len(female_chrs) == 2) & (chromosome in female_chrs):
        # ZW will trigger, XY won't
        return 1
    return 2


def get_input_vcf_for_genotype__variant_filtration(wildcards):
    region = wildcards.region
    return (
        POSTERIORS / f"{region}.vcf.gz"
        if region in DIPLOID_REGIONS
        else GENOTYPE_GVCFS / f"{region}.vcf.gz"
    )


def get_interval_for_haplotype_caller(wildcards):
    region = wildcards.region
    chrom, chrom_start, chrom_end, _ = REGIONS_BED4[REGIONS_BED4.name == region].values[
        0
    ]
    return f"{chrom}:{chrom_start}-{chrom_end}"
