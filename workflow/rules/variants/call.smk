include: "call_functions.smk"


rule variants__call__haplotype_caller:
    """Call variants for a single sample and region

    Note: if ploidy is correct, run it, if not (e.g. females have no Y chr) produce
    a file genotyping an empty interval.
    """
    input:
        reference=ancient(REFERENCE / f"{HOST_NAME}.fa.gz"),
        cram=RECALIBRATE / "{sample_id}.cram",
        crai=RECALIBRATE / "{sample_id}.cram.crai",
        dict_=REFERENCE / f"{HOST_NAME}.dict",
    output:
        gvcf_gz=temp(CALL / "{sample_id}" / "{region}.gvcf.gz"),
    log:
        CALL / "{sample_id}" / "{region}.log",
    benchmark:
        CALL / "{sample_id}" / "{region}.benchmark.tsv"
    retries: 5
    conda:
        "../../environments/gatk4.yml"
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=24 * 60,
    params:
        ploidy=get_ploidy_of_sample_and_chromosome,
        interval=get_interval_for_haplotype_caller,
        mock_interval=generate_mock_interval,
    shell:
        """
        if [[ {params.ploidy} -eq 0 ]]; then
            gatk HaplotypeCaller \
                --emit-ref-confidence GVCF \
                --input {input.cram} \
                --intervals {params.mock_interval} \
                --output {output.gvcf_gz} \
                --reference {input.reference} \
                --sample-ploidy 1 \
                2>{log} 1>&2
        else
            gatk HaplotypeCaller \
                --emit-ref-confidence GVCF \
                --input {input.cram} \
                --intervals {params.interval} \
                --output {output.gvcf_gz} \
                --reference {input.reference} \
                --sample-ploidy {params.ploidy} \
                2>{log} 1>&2
        fi
        """


rule variants__call__haplotype_caller__all:
    """Run HaplotypeCaller for all samples and regions"""
    input:
        [
            CALL / f"{sample_id}" / f"{region}.gvcf.gz"
            for sample_id in SAMPLES
            for region in REGIONS
        ],


rule variants__call__combine_gvcfs:
    """Combine gVCFs from multiple samples and one region"""
    input:
        gvcfs=get_files_to_genotype,
        ref=ancient(REFERENCE / f"{HOST_NAME}.fa.gz"),
    output:
        gvcf=temp(CALL / "{region}.vcf.gz"),
    log:
        CALL / "{region}.log",
    benchmark:
        CALL / "{region}.benchmark.tsv"
    resources:
        mem_mb=16 * 1024,
        runtime=24 * 60,
    wrapper:
        "v7.9.1/bio/gatk/combinegvcfs"


rule variants__call__combine_gvcfs__all:
    """Run CombineGVCFs for all regions"""
    input:
        [CALL / f"{region}.vcf.gz" for region in REGIONS],


rule variants__call__all:
    """Call variants for all libraries and chromosomes"""
    input:
        rules.variants__call__haplotype_caller__all.input,
        rules.variants__call__combine_gvcfs__all.input,
