use rule helpers__recompress as reference__recompress__fasta with:
    input:
        gz=ancient(features["dna"]),
    output:
        gz=REFERENCE / f"{HOST_NAME}.fa.gz",
    log:
        REFERENCE / f"{HOST_NAME}.fa.log",


use rule helpers__recompress as reference__recompress__vcf with:
    input:
        gz=ancient(features["known_vcf"]),
    output:
        gz=REFERENCE / f"{HOST_NAME}.vcf.gz",
    log:
        REFERENCE / f"{HOST_NAME}.vcf.log",


rule reference__recompress__gtf:
    """Sort and compress the GTF file from features.yaml with bedtools sort and bgzip"""
    input:
        gtf_gz=ancient(features["gtf"]),
    output:
        gtf_gz=REFERENCE / f"{HOST_NAME}.gtf.gz",
    log:
        REFERENCE / f"{HOST_NAME}.gtf.log",
    cache: "omit-software"
    conda:
        "../../environments/reference.yml"
    threads: 8
    resources:
        mem_mb=8 * 1024,
        runtime=1 * 60,
    shell:
        """
        (
            bedtools sort \
                -i {input.gtf_gz} \
                | bgzip \
                    --threads {threads} \
                    >{output.gtf_gz}
        ) 2>{log}
        """


rule reference__recompress__all:
    input:
        rules.reference__recompress__fasta.output,
        rules.reference__recompress__vcf.output,
        rules.reference__recompress__gtf.output,
