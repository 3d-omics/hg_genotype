def compose_rg_id(wildcards):
    """Compose read group ID for bowtie2"""
    return f"{wildcards.sample}_{wildcards.library}"


def compose_rg_extra(wildcards):
    """Compose read group extra information for bowtie2"""
    return f"LB:truseq_{wildcards.library}\tPL:Illumina\tSM:{wildcards.sample}"
