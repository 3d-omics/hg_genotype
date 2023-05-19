rule snpeff_download:
    """Download a SNPEff database"""
    output:
        db=directory("resources/snpeff/{snpeff_db}/"),
    log:
        "resources/snpeff/{snpeff_db}.log",
    params:
        datadir="$PWD/resources/snpeff/",
        snpeff_db="{snpeff_db}",
    conda:
        "../envs/snpeff.yml"
    shell:
        """
        snpEff download \
            {params.snpeff_db} \
            -dataDir {params.datadir} \
            -verbose \
        2> {log} 1>&2
        """


rule snpeff_ann:
    """Annotate variants with a SNPEff database"""
    input:
        vcf=GATK / "variants_filtered.vcf.gz",
        db="resources/snpeff/{snpeff_db}/",
    output:
        vcf=SNPEFF / "variants_{snpeff_db}.vcf.gz",
        genes=SNPEFF / "snpEff_stats_{snpeff_db}.genes.txt",
        html=SNPEFF / "snpEff_summary_{snpeff_db}.html",
        csv=SNPEFF / "snpEff_stats_{snpeff_db}.csv",
    log:
        SNPEFF / "variants_{snpeff_db}.log",
    conda:
        "../envs/snpeff.yml"
    params:
        snpeff_db="{snpeff_db}",
        datadir="$PWD/resources/snpeff/",
        html="snpEff_summary.html",
    resources:
        mem_mb=8000,
        runtime=24 * 60,
    shell:
        """
        (snpEff ann \
            {params.snpeff_db} \
            -dataDir {params.datadir} \
            -csvStats {output.csv} \
            -verbose \
            -i vcf \
            -o gatk \
            {input.vcf} \
        | bgzip \
            --compress-level 9 \
            --stdout \
        > {output.vcf} \
        ) 2> {log} 1>&2

        mv {params.html} {output.html} 2>> {log} 1>&2
        """


rule snpeff_report:
    """Generate SNPeff reports for all databases"""
    input:
        [SNPEFF / f"snpEff_stats_{db}.csv" for db in SNPEFF_DBS],


rule snpeff:
    """Run SNPeff for all databases"""
    input:
        rules.snpeff_report.input,
