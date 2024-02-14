include: "__functions__.smk"
include: "index.smk"
include: "map.smk"
include: "split.smk"
include: "mark_duplicates.smk"
include: "recalibrate.smk"
include: "merge.smk"


rule align:
    """Run all picard steps and get all reports"""
    input:
        rules.align__merge__all.input,
