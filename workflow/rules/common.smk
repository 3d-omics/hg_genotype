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
