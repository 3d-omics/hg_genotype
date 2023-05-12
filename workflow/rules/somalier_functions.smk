def compose_somalier_extract_one_param_out_dir(wildcards):
    """Compose the output directory for somalier extract"""
    return f"results/somalier/extract/{wildcards.sample}.{wildcards.library}"
