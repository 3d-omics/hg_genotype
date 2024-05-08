rule swaps__somalier__find_sites__:
    input:
        vcf=VARIANTS / "filter" / "all.filtered.vcf.gz",
        tbi=VARIANTS / "filter" / "all.filtered.vcf.gz.tbi",
    output:
        vcf=SOMALIER / "sites.vcf.gz",
    log:
        SOMALIER / "sites.log",
    conda:
        "__environment__.yml"
    params:
        min_allele_number=5,
        min_allele_frequency=0.15,
    shell:
        """
        somalier find-sites \
            --min-AN {params.min_allele_number} \
            --min-AF {params.min_allele_frequency} \
            --output-vcf {output.vcf} \
            {input.vcf} \
        2> {log} 1>&2
        """


rule swaps__somalier__extract__:
    input:
        sites=SOMALIER / "sites.vcf.gz",
        variants=VARIANTS / "filter" / "all.filtered.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
        fai=REFERENCE / "genome.fa.gz.fai",
    output:
        [SOMALIER / "extracted" / f"{sample_id}.somalier" for sample_id in SAMPLES],
    log:
        SOMALIER / "extracted.log",
    conda:
        "__environment__.yml"
    params:
        out_dir=SOMALIER / "extracted",
    shell:
        """
        somalier extract \
            --out-dir {params.out_dir} \
            --sites   {input.sites} \
            --fasta   {input.reference} \
            {input.variants} \
        2> {log} 1>&2
        """


rule swaps__somalier__relate__:
    input:
        extracted=[
            SOMALIER / "extracted" / f"{sample_id}.somalier" for sample_id in SAMPLES
        ],
    output:
        groups=SOMALIER / "relate.groups.tsv",
        html=SOMALIER / "relate.hmtl",
        pairs=SOMALIER / "relate.pairs.tsv",
        samples=SOMALIER / "relate.samples.tsv",
    log:
        SOMALIER / "relate.log",
    conda:
        "__environment__.log"
    params:
        output_prefix=SOMALIER / "relate",
    shell:
        """
        somalier relate \
            --output-prefix {params.output_prefix} \
            --infer \
            {input.extracted} \
        2> {log} 1>&2
        """


rule swaps__somalier__report:
    input:
        pairs=SOMALIER / "relate.pairs.tsv",
        samples=SOMALIER / "relate.samples.tsv",


rule swaps__somalier:
    input:
        rules.swaps__somalier__relate__.output,
