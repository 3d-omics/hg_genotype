rule gatk4_base_recalibrator:
    """
    TODO: add known sites
    """
    input:
        bam=PICARD / "{sample}.{library}.bam",
        bai=PICARD / "{sample}.{library}.bam.bai",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        known_sites=REFERENCE / "known_variants.vcf.gz",
        csi=REFERENCE / "known_variants.vcf.gz.tbi",
    output:
        table=GATK / "{sample}.{library}.base_recalibrator.txt",
    log:
        GATK / "{sample}.{library}.base_recalibrator.log",
    conda:
        "../envs/gatk4.yml"
    params:  # TODO: get from params.yaml
        extra="",
        java_options="",
    threads: 1
    shell:
        """
        gatk BaseRecalibrator \
            --input {input.bam} \
            --reference {input.reference} \
            --known-sites {input.known_sites} \
            --output {output.table} \
        2> {log} 1>&2
        """


rule gatk4_base_recalibrator_all:
    """Collect fasqtc reports from the results of fastp"""
    input:
        [
            GATK / f"{sample}.{library}.base_recalibrator.txt"
            for sample, library in SAMPLE_LIB
        ],


rule gatk4_apply_bqsr:
    input:
        bam=PICARD / "{sample}.{library}.bam",
        reference=REFERENCE / "genome.fa.gz",
        table=GATK / "{sample}.{library}.base_recalibrator.txt",
    output:
        bam=GATK / "{sample}.{library}.bqsr.bam",
    log:
        GATK / "{sample}.{library}.apply_bqsr.log",
    conda:
        "../envs/gatk4.yml"
    params:
        None,
    shell:
        """
        gatk ApplyBQSR \
            --input {input.bam} \
            --reference {input.reference} \
            --bqsr-recal-file {input.table} \
            --output {output.bam} \
        2> {log} 1>&2
        """


rule gatk4_apply_bqsr_all:
    """Collect fasqtc reports from the results of fastp"""
    input:
        [GATK / f"{sample}.{library}.bqsr.bam" for sample, library in SAMPLE_LIB],


# rule gatk4_analyze_covariates:  # TODO
#     shell:
#         pass


rule gatk4_haplotype_caller:  # TODO: parallelize this?
    input:
        reference=REFERENCE / "genome.fa.gz",
        bam=GATK / "{sample}.{library}.bqsr.bam",
    output:
        gvcf_gz=GATK / "{sample}.{library}.haplotype_caller.gvcf.gz",
    log:
        GATK / "{sample}.{library}.haplotype_caller.log",
    conda:
        "../envs/gatk4.yml"
    shell:
        """
        gatk HaplotypeCaller \
            --reference {input.reference} \
            --input {input.bam} \
            --output {output.gvcf_gz} \
            --emit-ref-confidence GVCF \
        2> {log} 1>&2
        """


rule gatk4_haplotype_caller_all:
    input:
        [
            GATK / f"{sample}.{library}.haplotype_caller.gvcf.gz"
            for sample, library in SAMPLE_LIB
        ],


# TODO: I think here comes the GenomicsDBImport step


rule gatk4_combine_gvcfs:
    input:
        vcf_gzs=get_files_to_genotype,
        reference=REFERENCE / "genome.fa.gz",
    output:
        vcf_gz=GATK / "joint_variants.vcf.gz",
    log:
        GATK / "joint_variants.log",
    conda:
        "../envs/gatk4.yml"
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


rule gatk4_genotype_gvcfs:
    input:
        vcf_gz=GATK / "joint_variants.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
    output:
        vcf_gz=GATK / "genotyped_variants.vcf.gz",
    log:
        GATK / "genotyped_variants.log",
    conda:
        "../envs/gatk4.yml"
    shell:
        """
        gatk GenotypeGVCFs \
            --variant {input.vcf_gz} \
            --reference {input.reference} \
            --output {output.vcf_gz} \
        2> {log} 1>&2
        """


# According to GATK Best Practices, we should be using VariantRecalibrator and
# Apply VQSR, but since they need datasets with known variants (human and
# genome-wide), we skip these steps.
# According to GATK, https://gatk.broadinstitute.org/hc/en-us/articles/360035532412-Can-t-use-VQSR-on-non-model-organism-or-small-dataset
# We have to :
#   - hard filter variants instead, and/or
#   - use CNNs to learn from the data
# So far, CNNs are only fit for single-sample SNPs, so we are voiding it


rule gatk4_calculate_genotype_posteriors:
    input:
        vcf=GATK / "genotyped_variants.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
    output:
        vcf=GATK / "variants_posteriors.vcf.gz",
    log:
        GATK / "variants_posteriors.log",
    conda:
        "../envs/gatk4.yml"
    shell:
        """
        gatk CalculateGenotypePosteriors \
            --output {output.vcf} \
            --variant {input.vcf} \
            --reference {input.reference} \
        2> {log} 1>&2
        """


rule gatk4_variant_filtration:
    input:
        reference=REFERENCE / "genome.fa.gz",
        vcf=GATK / "variants_posteriors.vcf.gz",
    output:
        vcf=GATK / "variants_filtered.vcf.gz",
    log:
        GATK / "variants_filtered.log",
    conda:
        "../envs/gatk4.yml"
    params:
        filter_name="HardFilter",
        filter_expression="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0",
    shell:
        """
        gatk VariantFiltration \
            --reference {input.reference} \
            --variant {input.vcf} \
            --output {output.vcf} \
            --filter-expression '{params.filter_expression}' \
            --filter-name '{params.filter_name}' \
        2> {log} 1>&2
        """


rule gatk4_variant_annotator:
    """
    Right now this step is useless
    TODO: use known_sites with --resource foo:resource.vcf
    TODO: use expressions with --expression foo.AF \
    TODO: add to params the annotations we want
    """
    input:
        bams=[GATK / f"{sample}.{library}.bqsr.bam" for sample, library in SAMPLE_LIB],
        bais=[
            GATK / f"{sample}.{library}.bqsr.bam.bai" for sample, library in SAMPLE_LIB
        ],
        vcf=GATK / "variants_filtered.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
        reference_vcf=REFERENCE / "known_variants.vcf.gz",
    output:
        vcf=GATK / "variants_annotated.vcf.gz",
    log:
        GATK / "variants_annotated.log",
    conda:
        "../envs/gatk4.yml"
    params:  # TODO: add annotations that we want
        input_bam=compose_input_bqsr_bams,
    shell:
        """
        gatk VariantAnnotator \
            --reference {input.reference} \
            {params.input_bam} \
            --variant {input.vcf} \
            --output {output.vcf} \
            --annotation Coverage \
        2> {log}
        """


# Now should come funcotator, but it is human centric
# We use snpeff instead
# rule gatk4_funcotator


# rule gatk_cnn_score_variants:
#     input:
#         bams = [
#             GATK / f"{sample}.{library}.bqsr.bam"
#             for sample, library in SAMPLE_LIB
#         ],
#         # TODO: add .bais here
#         reference = REFERENCE / "genome.fa.gz",
#         vcf =  GATK / "genotyped_variants.vcf.gz"
#     output:
#         vcf =  GATK / "annotated_variants.vcf.gz"
#     log:   GATK / "annotated_variants.log"
#     conda: "../envs/gatk4.yml"
#     params:
#         input_bams = compose_cnn_input_bams
#     shell:
#         """
#         gatk CNNScoreVariants \
#             {params.input_bams} \
#             --variant {input.vcf} \
#             --reference {input.reference} \
#             --output {output.vcf} \
#             --tensor-type read_tensor \
#         2> {log} 1>&2
#         """


rule gatk4_all:
    input:
        # rules.gatk4_base_recalibrator_all.input,
        # rules.gatk4_apply_bqsr_all.input,
        # rules.gatk4_haplotype_caller_all.input,
        # rules.gatk4_combine_gvcfs.output,
        # rules.gatk4_genotype_gvcfs.output,
        # GATK / "genotyped_variants.vcf.gz",
        # GATK / "variants_posteriors.vcf.gz",
        # GATK / "variants_filtered.vcf.gz",
        GATK / "variants_annotated.vcf.gz",
