BWAMEM2_INDEX_EXTENSIONS = ["amb", "bwt.2bit.64", "pac", "0123", "ann"]


def compose_read_group_header(wildcards):
    """Compose read group header for bwa mem"""

    rg_identifier = f"{wildcards.sample_id}_{wildcards.library_id}"
    rg_library = f"LB:truseq_{wildcards.library_id}"
    rg_platform = "PL:Illumina"
    rg_sample = f"SM:{wildcards.sample_id}"

    return f"@RG\\tID:{rg_identifier}\\t{rg_library}\\t{rg_platform}\\t{rg_sample}"
