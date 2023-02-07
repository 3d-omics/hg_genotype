rule fastp_sample_library:
    """Run fastp in default mode"""
    input:
        sample=[
            READS + "/{sample}.{library}_1.fq.gz",
            READS + "/{sample}.{library}_2.fq.gz",
        ],
    output:
        trimmed=[
            FASTP + "/{sample}.{library}_1.fq.gz",
            FASTP + "/{sample}.{library}_2.fq.gz",
        ],
        # Unpaired reads separately
        unpaired1=FASTP + "/{sample}.{library}_u1.fq.gz",
        unpaired2=FASTP + "/{sample}.{library}_u2.fq.gz",
        html=FASTP + "/{sample}.{library}.html",
        json=FASTP + "/{sample}.{library}.json",
    log:
        FASTP + "/{sample}.{library}.log",
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    threads: 2
    wrapper:
        "v1.23.1/bio/fastp"


rule fastp_trim:
    """Collect fastp files"""
    input:
        [
            f"{FASTP}/{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
        ],


rule fastp_fastqc:
    """Collect fasqtc reports from the results of fastp"""
    input:
        [
            f"{FASTP}/{sample}.{library}_{end}_fastqc.html"
            for sample, library in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
        ],


rule fastp:
    input:
        rules.fastp_trim.input,
        rules.fastp_fastqc.input,
