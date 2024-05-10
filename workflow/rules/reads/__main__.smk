include: "__functions__.smk"
include: "link.smk"
include: "fastqc.smk"


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads__link__all.input,
        rules.reads__fastqc__all.input,
