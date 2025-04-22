rule annotate__vep__tmp_vcf:
    """Slice the VCF file to only include the sample of interest."""
    input:
        FILTER / "all.filtered.vcf.gz",
    output:
        temp(VEP / "{sample}.bcf"),
    log:
        VEP / "{sample}.bcf.log",
    params:
        extra=lambda w: f"--samples {w.sample} --trim-alt-alleles",
    threads: 8
    wrapper:
        "v5.2.1/bio/bcftools/view"


rule annotate__vep__downlaod_plugins:
    """Download the VEP plugins"""
    output:
        directory(VEP / "plugins"),
    params:
        release=100,
    wrapper:
        "v5.2.1/bio/vep/plugins"


rule annotate__vep:
    """Annotate the VCF file with VEP"""
    input:
        calls=VEP / "{sample}.bcf",
        fasta=REFERENCE / f"{HOST_NAME}.fa.gz",
        gff=REFERENCE / f"{HOST_NAME}.gff.gz",
        gff_tbi=REFERENCE / f"{HOST_NAME}.gff.gz.tbi",
        plugins=VEP / "plugins",
    output:
        calls=VEP / "{sample}.vcf.gz",
        stats=VEP / "{sample}.vep.html",
    log:
        VEP / "{sample}.log",
    params:
        extra="--buffer_size 500",
        plugins=[],
    threads: 8
    resources:
        mem_mb=16 * 1024,
        runtime=8 * 60,
    wrapper:
        "v5.2.1/bio/vep/annotate"


rule annotate__vep__all:
    """Annotate all vcf files with VEP"""
    input:
        [VEP / f"{sample}.vcf.gz" for sample in SAMPLES],
