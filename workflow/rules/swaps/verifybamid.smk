rule swaps__verifybamid__generate__:
    input:
        vcf=REFERENCE / "known_variants.vcf.gz",
        fasta=REFERENCE / "genome.fa.gz",
    output:
        ud=SWAPS / "known_variants.UD",
        mu=SWAPS / "known_variants.mu",
        bed=SWAPS / "known_variants.bed",
        vcf=SWAPS / "known_variants.vcf.gz",
    params:
        out_prefix=SWAPS / "known_variants",
    log:
        SWAPS / "known_variants.log",
    conda:
        "__environment__.yml"
    shell:
        """
        rsync -Pravt {input.vcf} {output.vcf} 2> {log}

        verifybamid2 \
            --RefVCF    {output.vcf} \
            --Reference {input.fasta} \
        2>> {log} 1>&2
        """


rule swaps__verifybamid__run__:
    input:
        ud=SWAPS / "known_variants.UD",
        mu=SWAPS / "known_variants.mu",
        bed=SWAPS / "known_variants.bed",
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
