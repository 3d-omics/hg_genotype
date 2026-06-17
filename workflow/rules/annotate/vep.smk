rule annotate__vep__tmp_vcf:
    """Slice the VCF file to only include the sample of interest."""
    input:
        FILTER / "all.filtered.vcf.gz",
    output:
        temp(VEP / "{sample}.bcf"),
    log:
        VEP / "{sample}.bcf.log",
    threads: 2
    resources:
        mem_mb=1 * 1024,
        runtime=1 * 60,
    params:
        extra=lambda w: f"--samples {w.sample} --trim-alt-alleles",
    wrapper:
        "v7.9.1/bio/bcftools/view"


rule annotate__vep__download_cache:
    """Download the VEP cache for the reference species and build"""
    output:
        directory(VEP / "cache"),
    log:
        VEP / "cache.log",
    cache: "omit-software"
    params:
        species=features["species"],
        release=features["release"],
        build=HOST_NAME,
    wrapper:
        "v7.9.1/bio/vep/cache"


rule annotate__vep__download_plugins:
    """Download the VEP plugins"""
    output:
        directory(VEP / "plugins"),
    log:
        VEP / "plugins.log",
    cache: "omit-software"
    params:
        release=100,
    wrapper:
        "v7.9.1/bio/vep/plugins"


rule annotate__vep:
    """Annotate the VCF file with VEP"""
    input:
        calls=VEP / "{sample}.bcf",
        fasta=ancient(REFERENCE / f"{HOST_NAME}.fa.gz"),
        gtf=ancient(REFERENCE / f"{HOST_NAME}.gtf.gz"),
        gtf_tbi=REFERENCE / f"{HOST_NAME}.gtf.gz.tbi",
        cache=VEP / "cache",
        plugins=VEP / "plugins",
    output:
        calls=VEP / "{sample}.vcf.gz",
        stats=VEP / "{sample}.vep.html",
    log:
        VEP / "{sample}.log",
    resources:
        mem_mb=16 * 1024,
        runtime=8 * 60,
    params:
        extra="--buffer_size 500 --everything --warning_file /dev/stderr",
        plugins=[],
    wrapper:
        "v7.9.1/bio/vep/annotate"


rule annotate__vep__all:
    """Annotate all vcf files with VEP"""
    input:
        [VEP / f"{sample}.vcf.gz" for sample in SAMPLES],
