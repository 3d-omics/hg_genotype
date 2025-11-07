include: "mark_duplicates_functions.smk"


rule align__mark_duplicates:
    """Mark duplicates for all contigs and merging samples from different libraries"""
    input:
        bams=get_crams_for_mark_duplicates,
        crais=get_crais_for_mark_duplicates,
    output:
        bam=MARK_DUPLICATES / "{sample_id}.cram",
        metrics=MARK_DUPLICATES / "{sample_id}.metrics.tsv",
    log:
        MARK_DUPLICATES / "{sample_id}.bam.log",
    params:
        samtools_opts="--threads 24",
    threads: 24
    resources:
        mem_mb=8 * 1024,
        runtime=6 * 60,
    wrapper:
        "v5.2.1/bio/picard/markduplicates"
        # Do not update. samtools will output CRAM 3.1 and gatk/picard cannot
        # read it 2025-11-07


rule align__mark_duplicates__all:
    """Mark duplicates in all chromosomes and all libraries"""
    input:
        [MARK_DUPLICATES / f"{sample_id}.cram" for sample_id in SAMPLES],
