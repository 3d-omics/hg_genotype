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


def compose_rg_id(wildcards):
    return f"{wildcards.sample}_{wildcards.library}"


def compose_rg_extra(wildcards):
    return f"LB:truseq_{wildcards.library}\tPL:Illumina\tSM:{wildcards.sample}"


def get_files_to_genotype(wildcards):
    return [
        GATK / f"{sample}.{library}.haplotype_caller.gvcf.gz"
        for sample, library in SAMPLE_LIB
    ]


def compose_v_line(wildcards):
    files = [
        GATK / f"{sample}.{library}.haplotype_caller.gvcf.gz"
        for sample, library in SAMPLE_LIB
    ]
    text = ""
    for file in files:
        text += f"--variant {file} "
    return text


def compose_input_bqsr_bams(wildcards):
    bams = [GATK / f"{sample}.{library}.bqsr.bam" for sample, library in SAMPLE_LIB]
    result = "".join(f"--input {bam} " for bam in bams)
    return result
# def compose_cnn_input_bams(wildcards):
#     bams = [
#         GATK / f"{sample}.{library}.bqsr.bam"
#         for sample, library in SAMPLE_LIB
#     ]
#     result = "".join(f"--input {bam} " for bam in bams)
#     return result
