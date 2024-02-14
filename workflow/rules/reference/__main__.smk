include: "recompress.smk"


rule reference:
    """Re-bgzip the reference genome and known variants"""
    input:
        rules.reference_recompress_genome.output,
        rules.reference_recompress_vcf.output,
