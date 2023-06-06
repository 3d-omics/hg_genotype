rule picard_extract_one:
    """Extract a single chromosome from a CRAM file"""
    input:
        cram=BOWTIE2 / "{sample}.{library}.cram",
        crai=BOWTIE2 / "{sample}.{library}.cram.crai",
        reference=REFERENCE / "genome.fa.gz",
    output:
        bam=temp(PICARD / "extract/{sample}.{library}.{chromosome}.bam"),
    log:
        PICARD / "extract/{sample}.{library}.{chromosome}.log",
    conda:
        "../envs/samtools.yml"
    params:
        chromosome=lambda wildcards: f"{wildcards.chromosome}",
    resources:
        mem_mb=8000,
        runtime=360,
    shell:
        """
        samtools view \
            --bam \
            --uncompressed \
            --output {output.bam} \
            {input.cram} \
            {params.chromosome} \
        2> {log} 1>&2
        """


rule picard_extract_all:
    """Extract all chromosomes (the ones in features.yml) from all libraries file"""
    input:
        [
            PICARD / f"extract/{sample}.{library}/{chromosome}.bam"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
        ],


rule picard_markduplicates_one:
    """Mark duplicates in a single chromosome from a single library"""
    input:
        bam=PICARD / "extract/{sample}.{library}.{chromosome}.bam",
        reference=REFERENCE / "genome.fa.gz",
    output:
        bam=temp(PICARD / "markduplicates/{sample}.{library}.{chromosome}.bam"),
        metrics=PICARD / "markduplicates/{sample}.{library}.{chromosome}.metrics.tsv",
    log:
        PICARD / "markduplicates/{sample}.{library}.{chromosome}.log",
    benchmark:
        PICARD / "markduplicates/{sample}.{library}.{chromosome}.bmk"
    conda:
        "../envs/picard.yml"
    resources:
        mem_mb=8000,
        runtime=360,
    shell:
        """
        picard MarkDuplicates \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            --COMPRESSION_LEVEL 1 \
            --REFERENCE_SEQUENCE {input.reference} \
        2> {log} 1>&2
        """


rule picard_markduplicates_all:
    """Mark duplicates in all chromosomes and all libraries"""
    input:
        [
            PICARD / f"markduplicates/{sample}.{library}.{chromosome}.bam"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
        ],


rule picard_report_all:
    """Generate reports for all chromosomes and all libraries"""
    input:
        [
            PICARD / f"markduplicates/{sample}.{library}.{chromosome}.{report}"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
            for report in PICARD_REPORTS
        ],


rule picard:
    """Run all picard steps and get all reports"""
    input:
        rules.picard_markduplicates_all.input,
        rules.picard_report_all.input,
