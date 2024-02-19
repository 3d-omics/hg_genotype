rule reference__recompress__genome:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["dna"],
    output:
        fa_gz=REFERENCE / "genome.fa.gz",
    log:
        REFERENCE / "genome.log",
    conda:
        "__environment__.yml"
    threads: 8
    shell:
        """
        (gzip \
            --decompres \
            --stdout {input.fa_gz} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference__recompress__vcf:
    """Extract the vcf.gz on config.yaml into known_variants.vcf.gz with bgzip"""
    input:
        vcf_gz=features["reference"]["known_vcf"],
    output:
        vcf_gz=REFERENCE / "known_variants.vcf.gz",
    log:
        REFERENCE / "known_variants.log",
    conda:
        "__environment__.yml"
    threads: 8
    shell:
        """
        (gzip -dc {input.vcf_gz} \
        | bgzip --threads {threads} \
        > {output.vcf_gz}) \
        2> {log}
        """


rule reference__recompress:
    input:
        rules.reference__recompress__genome.output,
        rules.reference__recompress__vcf.output,
