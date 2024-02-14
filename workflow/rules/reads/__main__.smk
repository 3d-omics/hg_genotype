include: "__functions__.smk"
include: "link.smk"
include: "fastqc.smk"


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads_link_all.input,
        rules.reads_fastqc_all.input,


localrules:
    reads_link_,
