rule snpeff_download:
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
    input:
        vcf=GATK / "variants_annotated.vcf.gz",
        db="resources/snpeff/{snpeff_db}/",  # TODO: uncode this
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
    input:
        [SNPEFF / f"snpEff_stats_{db}.csv" for db in SNPEFF_DBS],


rule snpeff:
    input:
        [SNPEFF / f"variants_{db}.vcf.gz" for db in SNPEFF_DBS],
