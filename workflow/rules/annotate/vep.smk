rule annotate__vep__:
    input:
        vcf=VARIANT_FILTRATION / "variants_filtered.vcf.gz",
        fa=REFERENCE / "genome.fa.gz",
        gtf=REFERENCE / "annotation.gtf.gz",
        gtf_tbi=REFERENCE / "annotation.gtf.gz.tbi",
    output:
        vcf=VEP / "variants_annotated.vcf.gz",
    log:
        VEP / "variants_annotated.log",
    conda:
        "__environment__.yml"
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --fasta {input.fa} \
            --gtf {input.gtf} \
            --compress_output bgzip \
        2> {log} 1>&2
        """


rule annotate__vep:
    input:
        rules.annotate__vep__.output,
