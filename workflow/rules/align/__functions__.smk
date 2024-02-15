# def compose_rg_id(wildcards):
#     """Compose read group ID for bowtie2"""
#     return f"{wildcards.sample}_{wildcards.library}"


# def compose_rg_extra(wildcards):
#     """Compose read group extra information for bowtie2"""
#     return f"LB:truseq_{wildcards.library}\tPL:Illumina\tSM:{wildcards.sample}"


def compose_read_group_header(wildcards):
    """Compose read group header for bwa mem"""
    rg_identifier = f"{wildcards.sample}_{wildcards.library}"
    rg_library = f"LB:truseq_{wildcards.library}"
    rg_platform = "PL:Illumina"
    rg_sample = f"SM:{wildcards.sample}"
    return f"@RG\\tID:{rg_identifier}\\t{rg_library}\\t{rg_platform}\\t{rg_sample}"
