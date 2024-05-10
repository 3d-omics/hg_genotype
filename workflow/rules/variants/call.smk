rule variants__call__haplotype_caller__:
    """Call variants for a single library and chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        cram=RECALIBRATE / "{sample_id}.cram",
        crai=RECALIBRATE / "{sample_id}.cram.crai",
        dict_=REFERENCE / "genome.dict",
    output:
        gvcf_gz=CALL / "{sample_id}" / "{region}.gvcf.gz",
    log:
        CALL / "{sample_id}" / "{region}.log",
    conda:
        "__environment__.yml"
    params:
        ploidy=get_ploidy_of_sample_and_chromosome,
        interval=get_interval_for_haplotype_caller,
        mock_interval=generate_mock_interval,
    shell:
        """
        if [[ {params.ploidy} -eq 0 ]] ; then
            gatk HaplotypeCaller \
                --reference {input.reference} \
                --input {input.cram} \
                --output {output.gvcf_gz} \
                --emit-ref-confidence GVCF \
                --intervals {params.mock_interval} \
                --sample-ploidy 1 \
            2> {log} 1>&2
        else
            gatk HaplotypeCaller \
                --reference {input.reference} \
                --input {input.cram} \
                --output {output.gvcf_gz} \
                --emit-ref-confidence GVCF \
                --intervals {params.interval} \
                --sample-ploidy {params.ploidy} \
            2> {log} 1>&2
        fi
        """


rule variants__call__combine_gvcfs__:
    """Combine gVCFs from multiple samples and one region"""
    input:
        vcf_gzs=get_files_to_genotype,
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf_gz=CALL / "{region}.vcf.gz",
    log:
        CALL / "{region}.log",
    conda:
        "__environment__.yml"
    params:
        variant_line=compose_v_line,
    shell:
        """
        gatk CombineGVCFs \
            {params.variant_line} \
            --reference {input.reference} \
            --output {output.vcf_gz} \
        2> {log} 1>&2
        """


rule variants__call:
    """Call variants for all libraries and chromosomes"""
    input:
        [CALL / f"{region}.vcf.gz" for region in REGIONS],
