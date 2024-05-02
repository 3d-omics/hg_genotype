include: "snpeff.smk"
include: "vep.smk"


rule annotate:
    input:
        # rules.annotate__snpeff.input,
        rules.annotate__vep.input,
