include: "snpeff.smk"


rule annotate:
    input:
        rules.annotate__snpeff.input,
