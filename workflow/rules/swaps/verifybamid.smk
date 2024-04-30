rule swaps__verifybamid__generate__:
    input:
        vcf=VARIANT_FILTRATION / "variants_filtered.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
    output:
        ud=SWAPS / "variants_filtered.UD",
        mu=SWAPS / "variants_filtered.mu",
        bed=SWAPS / "variants_filtered.bed",
        vcf=SWAPS / "variants_filtered.vcf.gz",
    params:
        out_prefix=SWAPS / "variants_filtered",
    log:
        SWAPS / "variants_filtered.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ln --symbolic {input.vcf} {output.vcf}

        verifybamid2 \
            --RefVCF    {input.vcf} \
            --Reference {input.reference} \
        2> {log} 1>&2
        """


rule swaps__verifybamid__run__:
    input:
        ud=SWAPS / "variants_filtered.UD",
        mu=SWAPS / "variants_filtered.mu",
        bed=SWAPS / "variants_filtered.bed",
        cram=MAP / "{sample_id}.{library_id}.cram",
    output:
        SWAPS / "{sample_id}.{library_id}.selfSM",
    params:
        out_prefix=lambda w: SWAPS / f"{w.sample_id}.{w.library_id}",
    log:
        SWAPS / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        verifybamid2 \
            --UDFile   {input.ud} \
            --BamFile  {input.cram} \
            --MeanFile {input.mu} \
            --BEDFile  {input.bed} \
            --Output   {params.out_prefix} \
        2> {log} 1>&2
        """


rule swaps__verifybamid:
    input:
        [
            SWAPS / f"{sample_id}.{library_id}.selfSM"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
