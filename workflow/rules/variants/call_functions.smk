def get_value_for_sample_and_column(sample_id, column):
    """Look up the value in column `column` and sample `sample_id` samples.tsv."""

    if column not in samples.columns:
        raise ValueError(
            f"samples.tsv has no '{column}' column, but it is mandatory for "
            "every sample."
        )

    values = (
        samples.loc[samples["sample_id"] == sample_id, column]
        .dropna()
        .drop_duplicates()
        .values.tolist()
    )

    if not values or values[0] == "":
        raise ValueError(
            f"Sample '{sample_id}' has no value for the mandatory "
            f"'{column}' column in samples.tsv."
        )

    return int(values[0])


def get_sample_ploidy(sample_id):
    """Get the per-individual ploidy for a sample (usually 2)."""
    return get_value_for_sample_and_column(sample_id, "ploidy")


def get_pool_size(sample_id):
    """Get the number of pooled individuals for a sample (usually 1)."""
    return get_value_for_sample_and_column(sample_id, "pool_size")


def get_ploidy_of_sample_and_chromosome(wildcards):
    """Get the ploidy of a sample and chromosome.

    Every chromosome is treated as autosomal (ploidy * pool_size) unless it
    appears in features['reference']['organelles'] (→ pool_size) or in any
    sex's features['reference']['sex_ploidy'] table (→ pool_size *
    sex_ploidy[sex][chromosome], or 0 when the chromosome is absent for that
    sex, e.g. W in males).  Ploidy 0 causes the calling rule to emit a mock
    GVCF placeholder instead of running HaplotypeCaller for real.
    """
    sample_id = wildcards.sample_id
    chromosome = get_chromosome_from_region(wildcards.region)
    ploidy = get_sample_ploidy(sample_id)
    pool_size = get_pool_size(sample_id)

    if chromosome in features["reference"]["organelles"]:
        return pool_size

    sex_ploidy = features["reference"].get("sex_ploidy", {})
    known_sex_chroms = {c for sex_chroms in sex_ploidy.values() for c in sex_chroms}
    if chromosome in known_sex_chroms:
        sex = get_sex_from_sample(sample_id)
        if sex not in sex_ploidy:
            raise ValueError(
                f"Sample '{sample_id}' has sex '{sex}', which is not a key of "
                f"features['reference']['sex_ploidy'] ({list(sex_ploidy)}). "
                "Add it there."
            )
        return pool_size * sex_ploidy[sex].get(chromosome, 0)

    return ploidy * pool_size


def get_interval_for_haplotype_caller(wildcards):
    region = wildcards.region
    chrom, chrom_start, chrom_end, _ = REGIONS_BED4[REGIONS_BED4.name == region].values[
        0
    ]
    return f"{chrom}:{chrom_start}-{chrom_end}"


def generate_mock_interval(wildcards):
    """A trivial 1bp interval on the first autosomal region in the BED4,
    used as a placeholder for absent sex chromosomes (ploidy 0)."""
    organelles = features["reference"]["organelles"]
    sex_ploidy = features["reference"].get("sex_ploidy", {})
    known_sex_chroms = {c for sex_chroms in sex_ploidy.values() for c in sex_chroms}
    row = REGIONS_BED4[
        ~REGIONS_BED4.chrom.isin(organelles)
        & ~REGIONS_BED4.chrom.isin(known_sex_chroms)
    ].iloc[0]
    return f"{row.chrom}:{row.chromStart}-{row.chromEnd}"


def get_files_to_genotype(wildcards):
    """Get files to genotype for a sample and region"""
    return [
        CALL / sample_id / f"{wildcards.region}.gvcf.gz"
        for sample_id in SAMPLES
        # if wildcards.region in get_chromosomes_from_sample(sample_id)
    ]


def get_chromosome_from_region(region):
    """Get the chromosome from a region"""
    return REGIONS_BED4[REGIONS_BED4.name == region].chrom.values[0]


def get_sex_from_sample(sample):
    """From a sample name, get it's sex from the samples table"""
    sex = (
        samples[samples["sample_id"] == sample]["sex"]
        .drop_duplicates()
        .values.tolist()[0]
    )
    return sex
