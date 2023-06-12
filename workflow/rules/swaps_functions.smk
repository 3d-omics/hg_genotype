def get_bams_for_bcftools_call(wildcards):
    """Compose the list of BAM files to be used for bcftools call."""
    chromosome = wildcards.chromosome
    return [
        SWAPS / f"rename_library/{sample}.{library}.{chromosome}.bam"
        for sample, library in SAMPLE_LIB
    ]
