include: "reference/recompress.smk"


rule reference__all:
    """Re-bgzip the reference genome, known variants and gtf"""
    input:
        rules.reference__recompress__all.input,
