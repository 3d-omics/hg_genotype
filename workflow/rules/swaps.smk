rule swaps_rename_library_one:
    input:
        bam=PICARD / "markduplicates/{sample}.{library}.{chromosome}.bam",
        reference=REFERENCE / "genome.fa.gz",
    output:
        bam=SWAPS / "rename_library/{sample}.{library}.{chromosome}.bam",
    log:
        SWAPS / "rename_library/{sample}.{library}.{chromosome}.log",
    params:
        sample_library="{sample}.{library}",
    conda:
        "../envs/swaps.yml"
    shell:
        """
        picard AddOrReplaceReadGroups \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --RGLB {params.sample_library} \
            --RGPL illumina \
            --RGPU {params.sample_library} \
            --RGSM {params.sample_library} \
            --COMPRESSION_LEVEL 9 \
        2> {log} 1>&2
        """


rule swaps_rename_library_all:
    input:
        [
            SWAPS / f"rename_library/{sample}.{library}.{chromosome}.bam"
            for sample, library in SAMPLE_LIB
            for chromosome in get_sample_chromosomes(sample)
        ],


def get_bams_for_bcftools_call(wildcards):
    chromosome = wildcards.chromosome
    return [
        SWAPS / f"rename_library/{sample}.{library}.{chromosome}.bam"
        for sample, library in SAMPLE_LIB
    ]


rule swaps_call_one_chromosome:
    input:
        bams=get_bams_for_bcftools_call,
        reference=REFERENCE / "genome.fa.gz",
    output:
        vcf=SWAPS / "call/{chromosome}.vcf.gz",
    log:
        SWAPS / "call/{chromosome}.log",
    conda:
        "../envs/swaps.yml"
    shell:
        """
        (bcftools mpileup \
            --output-type z9 \
            --fasta-ref {input.reference} \
            {input} \
        | bcftools call \
            --variants-only \
            --multiallelic-caller \
            --output-type z9 \
            --output {output.vcf} ) \
        2> {log} 1>&2
        """


rule swaps_call_all:
    input:
        [SWAPS / f"bcftools_mpileup/{chromosome}.vcf.gz" for chromosome in CHROMOSOMES],


rule swaps_filter_one:
    input:
        vcf=SWAPS / "call/{chromosome}.vcf.gz",
    output:
        vcf=SWAPS / "filter/{chromosome}.vcf.gz",
    log:
        SWAPS / "filter/{chromosome}.log",
    conda:
        "../envs/swaps.yml"
    params:
        min_qual=30,
    shell:
        """
        bcftools filter \
            --include 'QUAL>{params.min_qual}' \
            --output-type z9 \
            --output {output.vcf} \
            {input.vcf} \
        2> {log} 1>&2
        """


rule swaps_filter_all:
    input:
        [SWAPS / f"call/{chromosome}.vcf.gz" for chromosome in CHROMOSOMES],


rule swaps_merge:
    input:
        vcf=[SWAPS / f"filter/{chromosome}.vcf.gz" for chromosome in CHROMOSOMES],
    output:
        vcf=SWAPS / "merge.vcf.gz",
    log:
        SWAPS / "merge.log",
    conda:
        "../envs/swaps.yml"
    shell:
        """
        bcftools concat \
            --output-type z9 \
            --output {output.vcf} \
            {input.vcf} \
        2> {log} 1>&2
        """


rule swaps_process:
    input:
        vcf=SWAPS / "merge.vcf.gz",
    output:
        tsv=SWAPS / "gtcheck.tsv",
    log:
        SWAPS / "gtcheck.log",
    conda:
        "../envs/swaps.yml"
    shell:
        """
        bcftools gtcheck \
            {input.vcf} \
        > {output.tsv} 2> {log}
        """


rule swaps_plot:
    input:
        tsv=SWAPS / "gtcheck.tsv",
    output:
        pdf=SWAPS / "gtcheck.pdf",
    log:
        SWAPS / "plot.log",
    conda:
        "../envs/swaps.yml"
    shell:
        """
        Rscript workflow/scripts/plot_gtcheck.R \
            --infile {input} \
            --outfile {output} \
        2> {log} 1>&2
        """


rule swaps:
    input:
        SWAPS / "gtcheck.pdf",
