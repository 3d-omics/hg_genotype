include: "__functions__.smk"
include: "sample.smk"
include: "step.smk"


rule report:
    input:
        rules.report__step.input,
        rules.report__sample.input,
