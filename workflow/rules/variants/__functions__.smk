# call
def get_sex_from_sample(sample):
    """From a sample name, get it's sex from the samples table"""
    sex = (
        samples[samples["sample_id"] == sample]["sex"]
        .drop_duplicates()
        .values.tolist()[0]
    )
    return sex


def get_chromosomes_from_sample(sample):
    """Get the chromosomes for a given sample"""
    sex = get_sex_from_sample(sample)
    autosomes = [str(x) for x in features["reference"]["autosomes"]]
    male_chromosomes = features["reference"]["male_chromosomes"]
    female_chromosomes = features["reference"]["female_chromosomes"]
    mitochondria = features["reference"]["mitochondria"]
    if sex == "male":
        return autosomes + male_chromosomes + mitochondria
    else:
        return autosomes + female_chromosomes + mitochondria


def get_chromosome_from_region(region):
    """Get the chromosome from a region"""
    return REGIONS_BED4[REGIONS_BED4.name == region].chrom.values[0]


def get_ploidy_of_sample_and_chromosome(wildcards):
    """Get the ploidy of a sample and chromosome"""
    sample_id = wildcards.sample_id
    region = wildcards.region
    chromosome = get_chromosome_from_region(region)
    sex = get_sex_from_sample(sample_id)

    autosomes = [str(x) for x in features["reference"]["autosomes"]]
    male_chromosomes = features["reference"]["male_chromosomes"]
    female_chromosomes = features["reference"]["female_chromosomes"]
    mitochondria = features["reference"]["mitochondria"]

    if chromosome in mitochondria:
        return 1
    if chromosome in autosomes:
        return 2
    if sex == "male" and len(male_chromosomes) == 2 and chromosome in male_chromosomes:
        return 1
    if sex == "male" and len(male_chromosomes) == 1 and chromosome in male_chromosomes:
        return 2
    if (
        sex == "female"
        and len(female_chromosomes) == 2
        and chromosome in female_chromosomes
    ):
        return 1
    if (
        sex == "female"
        and len(female_chromosomes) == 1
        and chromosome in female_chromosomes
    ):
        return 2
    return 0


def generate_mock_interval(wildcards):
    autosomes = [str(x) for x in features["reference"]["autosomes"]]
    male_chromosomes = features["reference"]["male_chromosomes"]
    female_chromosomes = features["reference"]["female_chromosomes"]
    mitochondria = features["reference"]["mitochondria"]
    chromosomes = autosomes + male_chromosomes + female_chromosomes + mitochondria
    return f"{chromosomes[0]}:1-1"


def get_interval_for_haplotype_caller(wildcards):
    region = wildcards.region
    chrom, chrom_start, chrom_end, _ = REGIONS_BED4[REGIONS_BED4.name == region].values[
        0
    ]
    return f"{chrom}:{chrom_start}-{chrom_end}"


def get_files_to_genotype(wildcards):
    """Get files to genotype for a sample, library and chromosome"""
    return [
        CALL / sample_id / f"{wildcards.region}.gvcf.gz"
        for sample_id in SAMPLES
        # if wildcards.region in get_chromosomes_from_sample(sample_id)
    ]


def compose_v_line(wildcards):
    """Compose the -v line for gatk4 genotype gvcfs"""
    files = get_files_to_genotype(wildcards)
    text = ""
    for file in files:
        text += f"--variant {file} "
    return text


# genotype
def compose_merge_vcfs_input_line(wildcards):
    files = [CALL / f"{region}.vcf.gz" for region in REGIONS]
    string = " ".join(f"--INPUT {file}" for file in files)
    return string


# filter
def get_input_vcf_for_genotype__variant_filtration(wildcards):
    region = wildcards.region
    return (
        POSTERIORS / f"{region}.vcf.gz"
        if region in DIPLOID_REGIONS
        else GENOTYPE / f"{region}.vcf.gz"
    )
