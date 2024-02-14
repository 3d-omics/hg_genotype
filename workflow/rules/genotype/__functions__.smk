def get_files_to_genotype(wildcards):
    """Get files to genotype for a sample, library and chromosome"""
    return [
        HAPLOTYPE_CALLER / f"{sample}.{library}" / f"{wildcards.chromosome}.gvcf.gz"
        for sample, library in SAMPLE_LIB
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
    sample = wildcards.sample
    library = wildcards.library
    chromosome = wildcards.chromosome
    sex = samples[(samples["sample"] == sample) & (samples["library"] == library)][
        ["sex"]
    ]

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
