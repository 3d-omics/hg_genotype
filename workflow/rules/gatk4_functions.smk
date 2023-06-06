def get_files_to_genotype(wildcards):
    """Get files to genotype for a sample, library and chromosome"""
    return [
        GATK / f"haplotype_caller/{sample}.{library}.{wildcards.chromosome}.gvcf.gz"
        for sample, library in SAMPLE_LIB
    ]


def compose_v_line(wildcards):
    """Compose the -v line for gatk4 genotype gvcfs"""
    files = [
        GATK / f"haplotype_caller/{sample}.{library}.{wildcards.chromosome}.gvcf.gz"
        for sample, library in SAMPLE_LIB
    ]
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

    male_chrs = features["male_chrs"]
    female_chrs = features["female_chrs"]
    mitochondria = features["mitochondria"]

    if chormosome in mitochondria:
        return 1
    if lenth(male_chrs) == 2 & chromosome in male_chrs:
        return 1
    if length(female_chrs) == 2 & chromsome in female_chrs:
        return 1
    return 2
