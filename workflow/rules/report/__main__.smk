include: "__functions__.smk"
include: "library.smk"
include: "step.smk"


rule report:
    input:
        rules.report__step.input,
        rules.report__library.input,
