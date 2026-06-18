_SEX_PLOIDY = features.get("sex_ploidy", {})
KNOWN_SEX_CHROMS = {c for chroms in _SEX_PLOIDY.values() for c in chroms}


def get_sample_ploidy(sample_id):
    """Get the per-individual ploidy for a sample (usually 2)."""
    return int(samples.loc[samples["sample_id"] == sample_id, "ploidy"].iloc[0])


def get_pool_size(sample_id):
    """Get the number of pooled individuals for a sample (usually 1)."""
    return int(samples.loc[samples["sample_id"] == sample_id, "pool_size"].iloc[0])


def get_ploidy_of_sample_and_chromosome(wildcards):
    """Get the ploidy of a sample and chromosome.

    Every chromosome is treated as autosomal (ploidy * pool_size) unless it
    appears in features['organelles'] (→ pool_size) or in any
    sex's features['sex_ploidy'] table (→ pool_size *
    sex_ploidy[sex][chromosome], or 0 when the chromosome is absent for that
    sex, e.g. W in males).  Ploidy 0 causes the calling rule to emit a mock
    GVCF placeholder instead of running HaplotypeCaller for real.
    """
    sample_id = wildcards.sample_id
    chromosome = get_chromosome_from_region(wildcards.region)
    ploidy = get_sample_ploidy(sample_id)
    pool_size = get_pool_size(sample_id)

    if chromosome in features["organelles"]:
        return pool_size

    if chromosome in KNOWN_SEX_CHROMS:
        sex = get_sex_from_sample(sample_id)
        if sex not in _SEX_PLOIDY:
            raise ValueError(
                f"Sample '{sample_id}' has sex '{sex}', which is not a key of "
                f"features['sex_ploidy'] ({list(_SEX_PLOIDY)}). "
                "Add it there."
            )
        return pool_size * _SEX_PLOIDY[sex].get(chromosome, 0)

    return ploidy * pool_size


def get_interval_for_haplotype_caller(wildcards):
    row = REGIONS_BED4[REGIONS_BED4.name == wildcards.region].iloc[0]
    return f"{row.chrom}:{row.chromStart}-{row.chromEnd}"


def generate_mock_interval(wildcards):
    """A trivial interval on the first autosomal region, used as placeholder for absent sex chromosomes."""
    row = REGIONS_BED4[
        ~REGIONS_BED4.chrom.isin(features["organelles"])
        & ~REGIONS_BED4.chrom.isin(KNOWN_SEX_CHROMS)
    ].iloc[0]
    return f"{row.chrom}:{row.chromStart}-{row.chromEnd}"


def get_files_to_genotype(wildcards):
    """Get files to genotype for a sample and region"""
    return [CALL / sample_id / f"{wildcards.region}.gvcf.gz" for sample_id in SAMPLES]


def get_chromosome_from_region(region):
    """Get the chromosome from a region."""
    return REGIONS_BED4[REGIONS_BED4.name == region]["chrom"].iloc[0]


def get_sex_from_sample(sample):
    """Get the sex of a sample from the samples table."""
    return samples[samples["sample_id"] == sample]["sex"].iloc[0]
