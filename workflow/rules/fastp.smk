rule fastp:
    """Run fastp in default mode"""
    input:
        forward_=READS + "{sample}.{library}_1.fq.gz",
        reverse_=READS + "{sample}.{library}_2.fq.gz",
    output:
        forward_=FASTP + "{sample}.{library}_1.fq.gz",
        reverse_=FASTP + "{sample}.{library}_2.fq.gz",
        unpaired1=FASTP + "{sample}.{library}_u1.fq.gz",
        unpaired2=FASTP + "{sample}.{library}_u2.fq.gz",
        html=FASTP + "{sample}.{library}.html",
        json=FASTP + "{sample}.{library}.json",
    log:
        FASTP + "{sample}.{library}.log",
    params:
        adapter_forward="ACGGCTAGCTA",  # TODO: get from samples.tsv
        adapter_reverse="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",  # TODO: get from samples.tsv
        extra_args="",  # TODO: fill with arguments from params.yml
    threads: MAX_THREADS
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 {output.forward_} \
            --out2 {output.reverse_} \
            --unpaired1 {output.unpaired1} \
            --unpaired2 {output.unpaired2} \
            --html {output.html} \
            --json {output.json} \
            --compression 1 \
            --verbose \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            {params.extra_args} \
            --thread {threads} \
        2> {log} 1>&2
        """


rule fastp_all_samples:
    """Collect fastp files"""
    input:
        [
            FASTP + f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
        ],


rule fastp_fastqc:
    """Collect fasqtc reports from the results of fastp"""
    input:
        [
            FASTP + f"{sample}.{library}_{end}_fastqc.html"
            for sample, library in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
        ],


rule fastp_all:
    input:
        rules.fastp_all_samples.input,
        rules.fastp_fastqc.input,
