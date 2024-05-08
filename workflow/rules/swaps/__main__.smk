include: "somalier.smk"
include: "verifybamid.smk"


rule swaps:
    input:
        rules.swaps__somalier.input,
