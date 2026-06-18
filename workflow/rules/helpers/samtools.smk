rule helpers__samtools__stats_cram:
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
        reference=ancient(REFERENCE / f"{HOST_NAME}.fa.gz"),
        fai=REFERENCE / f"{HOST_NAME}.fa.gz.fai",
        gzi=REFERENCE / f"{HOST_NAME}.fa.gz.gzi",
    output:
        "{prefix}.stats",
    log:
        "{prefix}.stats.log",
    conda:
        "../../environments/samtools.yml"
    resources:
        mem_mb=1 * 1024,
        runtime=1 * 60,
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
            >{output} \
            2>{log}
        """
