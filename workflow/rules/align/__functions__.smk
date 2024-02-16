def compose_read_group_header(wildcards):
    """Compose read group header for bwa mem"""
    rg_identifier = f"{wildcards.sample}_{wildcards.library}"
    rg_library = f"LB:truseq_{wildcards.library}"
    rg_platform = "PL:Illumina"
    rg_sample = f"SM:{wildcards.sample}"
    return f"@RG\\tID:{rg_identifier}\\t{rg_library}\\t{rg_platform}\\t{rg_sample}"


def get_crams_for_mark_duplicates(wildcards):
    """Get all the bams for mark duplicates. They will be joined by sample"""
    libraries = samples[samples["sample"] == wildcards.sample].library.values.tolist()
    crams = [MAP / f"{wildcards.sample}.{library}.cram" for library in libraries]
    return crams


def compose_input_line_for_mark_duplicates(wildcards):
    """Compose input line for mark duplicates"""
    crams = get_crams_for_mark_duplicates(wildcards)
    return " ".join(f"--INPUT {cram}" for cram in crams)
