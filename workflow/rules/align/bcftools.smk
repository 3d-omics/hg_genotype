include: "bcftools_functions.smk"


rule align__bcftools__call:
    input:
        crams=[MARK_DUPLICATES / f"{sample_id}.cram" for sample_id in SAMPLES],
        crais=[MARK_DUPLICATES / f"{sample_id}.cram.crai" for sample_id in SAMPLES],
        fasta=ancient(REFERENCE / f"{HOST_NAME}.fa.gz"),
        fai=REFERENCE / f"{HOST_NAME}.fa.gz.fai",
    output:
        bcf=temp(BCFTOOLS / "{window}.bcf"),
    log:
        BCFTOOLS / "{window}.log",
    conda:
        "../../environments/bcftools.yml"
    resources:
        mem_mb=2 * 1024,
        runtime=2 * 60,
    params:
        padded=get_bcftools_window_padded,
        core=get_bcftools_window_core,
    shell:
        """
        (
            bcftools mpileup \
                --fasta-ref {input.fasta} \
                --region {params.padded} \
                --output-type u \
                {input.crams} \
                | bcftools call \
                    --multiallelic-caller \
                    --variants-only \
                    --output-type u \
                | bcftools filter \
                    --include 'QUAL > 30' \
                    --targets {params.core} \
                    --targets-overlap 0 \
                    --output-type b \
                    --output {output.bcf}
        ) 2>{log}
        """


rule align__bcftools__concat:
    input:
        [BCFTOOLS / f"{window}.bcf" for window in BCFTOOLS_WINDOW_NAMES],
    output:
        BCFTOOLS / "bcftools.vcf.gz",
    log:
        BCFTOOLS / "bcftools.log",
    conda:
        "../../environments/bcftools.yml"
    resources:
        mem_mb=1 * 1024,
        runtime=1 * 60,
    shell:
        """
        bcftools concat \
            --output-type z \
            --output {output} \
            {input} \
            2>{log}
        """


rule align__bcftools__all:
    input:
        BCFTOOLS / "bcftools.vcf.gz",
