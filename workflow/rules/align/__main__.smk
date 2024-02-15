include: "__functions__.smk"
include: "index.smk"
include: "map.smk"
include: "mark_duplicates.smk"
include: "recalibrate.smk"


rule align:
    """Run all picard steps and get all reports"""
    input:
        rules.align__recalibrate__applybqsr__all.input,
