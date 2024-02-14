rule align__split__:
    """Extract a single chromosome from a CRAM file"""
    input:
        cram=MAP / "{sample}.{library}.cram",
        crai=MAP / "{sample}.{library}.cram.crai",
        reference=REFERENCE / "genome.fa.gz",
    output:
        bam=temp(SPLIT / "{sample}.{library}" / "{chromosome}.bam"),
    log:
        SPLIT / "{sample}.{library}" / "{chromosome}.log",
    conda:
        "__environment__.yml"
    params:
        chromosome=lambda w: f"{w.chromosome}",
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


rule align__split__all:
    """Extract all chromosomes (the ones in features.yml) from all libraries file"""
    input:
        [
            SPLIT / f"{sample}.{library}" / f"{chromosome}.bam"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
        ],
