def get_reads(wildcards):
    forward_, reverse_ = samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ][["forward", "reverse"]].values[0]
    return forward_, reverse_


def get_forward(wildcards):
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward"].tolist()[0]


def get_reverse(wildcards):
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse"].tolist()[0]


def get_forward_adapter(wildcards):
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward_adapter"].tolist()[0]


def get_reverse_adapter(wildcards):
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse_adapter"].tolist()[0]


def compose_rg_id(wildcards):
    return f"{wildcards.sample}_{wildcards.library}"


def compose_rg_extra(wildcards):
    return f"LB:truseq_{wildcards.library}\tPL:Illumina\tSM:{wildcards.sample}"


def get_files_to_genotype(wildcards):
    return [
        GATK / f"haplotype_caller/{sample}.{library}.{wildcards.chromosome}.gvcf.gz"
        for sample, library in SAMPLE_LIB
    ]


def compose_v_line(wildcards):
    files = [
        GATK / f"haplotype_caller/{sample}.{library}.{wildcards.chromosome}.gvcf.gz"
        for sample, library in SAMPLE_LIB
    ]
    text = ""
    for file in files:
        text += f"--variant {file} "
    return text


def get_picard_per_sample_files(wildcards):
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
    files = [
        GATK
        / f"base_recalibrator/{wildcards.sample}.{wildcards.library}.{chromosome}.txt"
        for chromosome in CHROMOSOMES
    ]
    return files


def compose_somalier_extract_one_param_out_dir(wildcards):
    return f"results/somalier/extract/{wildcards.sample}.{wildcards.library}"
