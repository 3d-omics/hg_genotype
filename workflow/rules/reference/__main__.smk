include: "recompress.smk"


rule reference:
    """Re-bgzip the reference genome and known variants"""
    input:
        rules.reference__recompress.input,
