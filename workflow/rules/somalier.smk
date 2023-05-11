rule somalier_extract_one:
    """Extract SNPs from one sample"""
    input:
        cram=BOWTIE2 / "{sample}.{library}.cram",
        crai=BOWTIE2 / "{sample}.{library}.cram.crai",
        reference=REFERENCE / "genome.fa.gz",
        vcf=REFERENCE / "known_variants.vcf.gz",
    output:
        somalier=SOMALIER / "extract/{sample}.{library}/{sample}.somalier",
    log:
        SOMALIER / "extract/{sample}.{library}.log",
    conda:
        "../envs/somalier.yml"
    params:
        out_dir=compose_somalier_extract_one_param_out_dir,
        sample_prefix="{sample}.{library}",  # libpath + wildcards don't work
    shell:
        """
        somalier extract \
            --sites {input.vcf} \
            --fasta {input.reference} \
            --out-dir {params.out_dir} \
            --sample-prefix {params.sample_prefix} \
            {input.cram} \
        2> {log} 1>&2
        """


rule somalier_extract_all:
    """Extract SNPs from all samples"""
    input:
        [
            SOMALIER / f"extract/{sample}.{library}/{sample}.somalier"
            for sample, library in SAMPLE_LIB
        ],


rule somalier_relate:
    """Relate all samples"""
    input:
        somalier=rules.somalier_extract_all.input,
    output:
        html=SOMALIER / "relate.html",
        groups=SOMALIER / "relate.groups.tsv",
        pairs=SOMALIER / "relate.pairs.tsv",
        samples=SOMALIER / "relate.samples.tsv",
    log:
        SOMALIER / "somalier.log",
    conda:
        "../envs/somalier.yml"
    params:
        min_depth=7,
        output_prefix=SOMALIER / "relate",
    shell:
        """
        somalier relate \
            --min-depth {params.min_depth} \
            --infer \
            --output-prefix {params.output_prefix} \
            {input.somalier} \
        2> {log} 1>&2
        """


rule somalier_report:
    """Generate report"""
    input:
        rules.somalier_relate.output,


rule somalier:
    """Run somalier"""
    input:
        rules.somalier_relate.output,
