SNPEFF_DB = features["snpeff"]["genome"]


rule snpeff_download:
    output:
        db=directory("resources/snpeff/{snpeff_db}/"),
        # TODO: ask for concrete files
        # TODO: un-hardcode database
    log:
        "resources/snpeff/snpEff_v5_1{snpeff_db}.log",
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
        vcf=GATK + "variants_annotated.vcf.gz",
        db="resources/snpeff/{snpeff_db}/",  # TODO: uncode this
    output:
        vcf=SNPEFF + "variants_{snpeff_db}.vcf.gz",
        genes=SNPEFF + "snpEff_genes_{snpeff_db}.txt",
        html=SNPEFF + "snpEff_summary_{snpeff_db}.html",
    log:
        SNPEFF + "variants_{snpeff_db}.log",
    conda:
        "../envs/snpeff.yml"
    params:
        snpeff_db="{snpeff_db}",
        datadir="$PWD/resources/snpeff/",
        genes="snpEff_genes.txt",
        html="snpEff_summary.html",
    shell:
        """
        (snpEff ann \
            {params.snpeff_db} \
            -dataDir {params.datadir} \
            -verbose \
            -i vcf \
            -o gatk \
            {input.vcf} \
        | bgzip \
            --compress-level 9 \
            --stdout \
        > {output.vcf} \
        ) 2> {log} 1>&2

        mv {params.genes} {output.genes} 2>> {log} 1>&2
        mv {params.html} {output.html} 2>> {log} 1>&2
        """


rule snpeff:
    input:
        SNPEFF + f"variants_{SNPEFF_DB}.vcf.gz",
