rule bowtie2_build:
    input:
        reference=REFERENCE / "genome.fa.gz",
    output:
        multiext(
            f"{REFERENCE}/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        BOWTIE2 / "build.log",
    conda:
        "../envs/bowtie2.yml"
    params:
        output_path=REFERENCE / "genome",
        extra="",  # optional parameters
    threads: 8
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.reference} \
            {params.output_path} \
        2> {log} 1>&2
        """


rule bowtie2_map:
    """Not using wrapper to sort it as it comes out"""
    input:
        forward_=FASTP / "{sample}.{library}_1.fq.gz",
        reverse_=FASTP / "{sample}.{library}_2.fq.gz",
        unpaired1=FASTP / "{sample}.{library}_u1.fq.gz",
        unpaired2=FASTP / "{sample}.{library}_u2.fq.gz",
        idx=multiext(
            f"{REFERENCE}/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        reference=REFERENCE / "genome.fa.gz",
    output:
        cram=BOWTIE2 / "{sample}.{library}.cram",
    log:
        BOWTIE2 / "{sample}.{library}.log",
    params:
        index_prefix=REFERENCE / "genome",
        extra="",  # optional parameters
        samtools_mem_per_thread="1G",
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: MAX_THREADS
    conda:
        "../envs/bowtie2.yml"
    shell:
        """
        (bowtie2 \
            -x {params.index_prefix} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            -U {input.unpaired1},{input.unpaired2} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -m {params.samtools_mem_per_thread} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2:
    input:
        [BOWTIE2 / f"{sample}.{library}.cram" for sample, library in SAMPLE_LIB],


rule bowtie2_reports:
    input:
        [
            BOWTIE2 / f"{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in "stats.tsv flagstats.txt idxstats.tsv".split()
        ],
