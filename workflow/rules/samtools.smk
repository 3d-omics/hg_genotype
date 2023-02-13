rule bai:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "{prefix}.bam.bai.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule crai:
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    log:
        "{prefix}.cram.crai.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule dict_fa:
    input:
        "{prefix}.fa",
    output:
        "{prefix}.dict",
    log:
        "{prefix}.dict.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools dict {input} --output {output} 2> {log} 1>&2"


rule dict_fagz:
    input:
        "{prefix}.fa.gz",
    output:
        "{prefix}.dict",
    log:
        "{prefix}.dict.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools dict {input} --output {output} 2> {log} 1>&2"


rule vcf_gz_tbi:
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    log:
        "{prefix}.vcf.gz.tbi.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "tabix {input}"


rule vcf_gz:
    input:
        "{prefix}.vcf",
    output:
        "{prefix}.vcf.gz",
    log:
        "{prefix}.vcf.gz.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "bgzip {input}"
