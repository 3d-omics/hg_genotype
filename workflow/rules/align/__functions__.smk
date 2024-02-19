def compose_read_group_header(wildcards):
    """Compose read group header for bwa mem"""
    rg_identifier = f"{wildcards.sample_id}_{wildcards.library_id}"
    rg_library = f"LB:truseq_{wildcards.library_id}"
    rg_platform = "PL:Illumina"
    rg_sample = f"SM:{wildcards.sample_id}"
    return f"@RG\\tID:{rg_identifier}\\t{rg_library}\\t{rg_platform}\\t{rg_sample}"


def get_crams_for_mark_duplicates(wildcards):
    """Get all the bams for mark duplicates. They will be joined by sample"""
    libraries = samples[
        samples["sample_id"] == wildcards.sample_id
    ].library_id.values.tolist()
    crams = [
        MAP / f"{wildcards.sample_id}.{library_id}.cram" for library_id in libraries
    ]
    return crams


def compose_input_line_for_mark_duplicates(wildcards):
    """Compose input line for mark duplicates"""
    crams = get_crams_for_mark_duplicates(wildcards)
    return " ".join(f"--INPUT {cram}" for cram in crams)
