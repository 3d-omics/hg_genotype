def get_reads(wildcards):
    """Get reads for a sample and library."""
    forward_, reverse_ = samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ][["forward", "reverse"]].values[0]
    return forward_, reverse_


def get_forward(wildcards):
    """Get forward reads for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward"].tolist()[0]


def get_reverse(wildcards):
    """Get reverse reads for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse"].tolist()[0]


def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward_adapter"].tolist()[0]


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse_adapter"].tolist()[0]


def compose_rg_id(wildcards):
    """Compose read group ID for bowtie2"""
    return f"{wildcards.sample}_{wildcards.library}"


def compose_rg_extra(wildcards):
    """Compose read group extra information for bowtie2"""
    return f"LB:truseq_{wildcards.library}\tPL:Illumina\tSM:{wildcards.sample}"


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


def get_picard_per_sample_files(wildcards):
    """Get picard reports for all samples, libraries, and chromosomes"""
    sample = wildcards.sample
    library = wildcards.library
    ANALYSES = ["stats.tsv", "flagstats.txt", "idxstats.tsv", "metrics.tsv"]
    files = [
        PICARD / f"markduplicates/{sample}.{library}.{chromosome}.{analysis}"
        for chromosome in CHROMOSOMES
        for analysis in ANALYSES
    ]
    return files


def get_gatk4_base_recalibrator_per_sample_files(wildcards):
    """Get gatk4 base recalibrator reports for all samples, libraries, and chromosomes"""
    files = [
        GATK
        / f"base_recalibrator/{wildcards.sample}.{wildcards.library}.{chromosome}.txt"
        for chromosome in CHROMOSOMES
    ]
    return files


def compose_somalier_extract_one_param_out_dir(wildcards):
    """Compose the output directory for somalier extract"""
    return f"results/somalier/extract/{wildcards.sample}.{wildcards.library}"
