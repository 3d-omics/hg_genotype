include: "__functions__.smk"
include: "call.smk"
include: "genotype.smk"
include: "filter.smk"


rule variants:
    input:
        rules.variants__filter.input,
