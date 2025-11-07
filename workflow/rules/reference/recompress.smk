rule reference__recompress__genome:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["dna"],
    output:
        fa_gz=REFERENCE / f"{HOST_NAME}.fa.gz",
    log:
        REFERENCE / f"{HOST_NAME}.fa.log",
    conda:
        "../../environments/reference.yml"
    cache: "omit-software"
    threads: 8
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.fa_gz} \
        | bgzip \
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
        vcf_gz=REFERENCE / f"{HOST_NAME}.vcf.gz",
    log:
        REFERENCE / f"{HOST_NAME}.vcf.log",
    conda:
        "../../environments/reference.yml"
    cache: "omit-software"
    threads: 8
    shell:
        """
        ( gzip \
            --decompress \
            --stdout \
            {input.vcf_gz} \
        | bgzip \
            --threads {threads} \
        > {output.vcf_gz}) \
        2> {log}
        """


rule reference__recompress__gff:
    """Sort and compress the GFF file from config.yaml into a .gff.gz file with bedtools sort and bgzip"""
    input:
        gff_gz=features["reference"]["gff"],
    output:
        gff_gz=REFERENCE / f"{HOST_NAME}.gff.gz",
    log:
        REFERENCE / f"{HOST_NAME}.gff.log",
    conda:
        "../../environments/reference.yml"
    cache: "omit-software"
    threads: 8
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        ( bedtools sort \
            -i {input.gff_gz} \
        | bgzip \
            --threads {threads} \
        > {output.gff_gz} \
        ) 2> {log}
        """


rule reference__recompress__all:
    input:
        rules.reference__recompress__genome.output,
        rules.reference__recompress__vcf.output,
        rules.reference__recompress__gff.output,
