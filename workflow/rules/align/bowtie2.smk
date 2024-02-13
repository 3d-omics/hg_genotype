rule bowtie2_build:
    """Build bowtie2 index

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
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
    benchmark:
        BOWTIE2 / "build.bmk"
    conda:
        "__environment__.yml"
    params:
        output_path=REFERENCE / "genome",
        extra=params["bowtie2"]["extra"],
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


rule bowtie2_map_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
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
        cram=protected(BOWTIE2 / "{sample}.{library}.cram"),
    log:
        BOWTIE2 / "{sample}.{library}.log",
    benchmark:
        BOWTIE2 / "{sample}.{library}.bmk"
    params:
        index_prefix=REFERENCE / "genome",
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
    conda:
        "__environment__.yml"
    resources:
        mem_mb=30000,
        runtime=1440,
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
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2_map_all:
    """Collect the results of `bowtie2_map_one` for all libraries"""
    input:
        [BOWTIE2 / f"{sample}.{library}.cram" for sample, library in SAMPLE_LIB],


rule bowtie2_report_all:
    """Generate bowtie2 reports for all libraries:
    - samtools stats
    - samtools flagstats
    - samtools idxstats
    """
    input:
        [
            BOWTIE2 / f"{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
        ],


rule bowtie2:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.bowtie2_map_all.input,
        rules.bowtie2_report_all.input,
