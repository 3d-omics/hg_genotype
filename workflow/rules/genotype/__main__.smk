include: "__functions__.smk"
include: "haplotype_caller.smk"
include: "combine_gvcfs.smk"
include: "genotype_gvcfs.smk"
include: "calculate_genotype_posteriors.smk"
include: "variant_filtration.smk"


rule genotype:
    input:
        rules.genotype__variant_filtration__all.input,
