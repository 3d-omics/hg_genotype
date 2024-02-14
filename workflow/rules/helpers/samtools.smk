rule helpers__samtools__bai__:
    """Generate a bam index"""
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "{prefix}.bam.bai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule helpers__samtools__crai__:
    """Generate a cram index"""
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    log:
        "{prefix}.cram.crai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule helpers__samtools__dict_fa__:
    """Generate a dictionary from a .fa"""
    input:
        "{prefix}.fa",
    output:
        "{prefix}.dict",
    log:
        "{prefix}.dict.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools dict {input} --output {output} 2> {log} 1>&2"


rule helpers__samtools__dict_fagz__:
    """Generate a dictionary from a fa.gz"""
    input:
        "{prefix}.fa.gz",
    output:
        "{prefix}.dict",
    log:
        "{prefix}.dict.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools dict {input} --output {output} 2> {log} 1>&2"


rule helpers__samtools__vcf_gz_tbi__:
    """Generate index for a vcf.gz"""
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    log:
        "{prefix}.vcf.gz.tbi.log",
    conda:
        "__environment__.yml"
    shell:
        "tabix {input} 2> {log} 1>&2"


rule helpers__samtools__vcf_gz__:
    """bgzip a vcf file"""
    input:
        "{prefix}.vcf",
    output:
        "{prefix}.vcf.gz",
    log:
        "{prefix}.vcf.gz.log",
    conda:
        "__environment__.yml"
    shell:
        "bgzip {input} 2> {log} 1>&2"


rule helpers__samtools__samtools_stats_bam__:
    """Compute stats for a bam"""
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai",
    output:
        tsv="{prefix}.stats.tsv",
    log:
        "{prefix}.stats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools stats {input.bam} > {output.tsv} 2> {log}"


rule helpers__samtools__samtools_stats_cram__:
    """Compute stats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.stats.tsv",
    log:
        "{prefix}.stats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools stats {input.cram} > {output.tsv} 2> {log}"


rule helpers__samtools__samtools_flagstats_bam__:
    """Compute flagstats for a bam"""
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools flagstats {input.bam} > {output.txt} 2> {log}"


rule helpers__samtools__samtools_flagstats_cram__:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule helpers__samtools__samtools_idxstats_bam__:
    """Compute idxstats for a bam"""
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools idxstats {input.bam} > {output.tsv} 2> {log}"


rule helpers__samtools__samtools_idxstats_cram__:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"
