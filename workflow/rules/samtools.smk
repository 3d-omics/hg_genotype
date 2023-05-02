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
        "tabix {input} 2> {log} 1>&2"


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
        "bgzip {input} 2> {log} 1>&2"


rule samtools_stats_bam:
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai",
    output:
        tsv="{prefix}.stats.tsv",
    log:
        "{prefix}.stats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools stats {input.bam} > {output.tsv} 2> {log}"


rule samtools_stats_cram:
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.stats.tsv",
    log:
        "{prefix}.stats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools stats {input.cram} > {output.tsv} 2> {log}"


rule samtools_flagstats_bam:
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools flagstats {input.bam} > {output.txt} 2> {log}"


rule samtools_flagstats_cram:
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule samtools_idxstats_bam:
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools idxstats {input.bam} > {output.tsv} 2> {log}"


rule samtools_idxstats_cram:
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"
