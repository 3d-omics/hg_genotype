rule gatk4_base_recalibrator_one:
    """Compute the recalibration table for a single library and chromosome"""
    input:
        bam=PICARD / "markduplicates/{sample}.{library}.{chromosome}.bam",
        bai=PICARD / "markduplicates/{sample}.{library}.{chromosome}.bam.bai",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        known_sites=REFERENCE / "known_variants.vcf.gz",
        csi=REFERENCE / "known_variants.vcf.gz.tbi",
    output:
        table=GATK / "base_recalibrator/{sample}.{library}.{chromosome}.txt",
    log:
        GATK / "base_recalibrator/{sample}.{library}.{chromosome}.log",
    benchmark:
        GATK / "base_recalibrator/{sample}.{library}.{chromosome}.bmk"
    conda:
        "../envs/gatk4.yml"
    params:
        extra=params["gatk4"]["base_recalibrator"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk BaseRecalibrator \
            {params.extra} \
            --input {input.bam} \
            --reference {input.reference} \
            --known-sites {input.known_sites} \
            --output {output.table} \
        2> {log} 1>&2
        """


rule gatk4_base_recalibrator_all:
    """Compute recalibration for all chromosomes and libraries"""
    input:
        [
            GATK / f"base_recalibrator/{sample}.{library}.{chromosome}.txt"
            for sample, library in SAMPLE_LIB
            for chromosome in CHROMOSOMES
        ],


rule gatk4_apply_bqsr_one:
    """Apply the recalibration table to a single library and chromosome"""
    input:
        bam=PICARD / "markduplicates/{sample}.{library}.{chromosome}.bam",
        reference=REFERENCE / "genome.fa.gz",
        table=GATK / "base_recalibrator/{sample}.{library}.{chromosome}.txt",
        dict_=REFERENCE / "genome.dict",
    output:
        bam=protected(GATK / "apply_bqsr/{sample}.{library}.{chromosome}.bam"),
    log:
        GATK / "apply_bqsr/{sample}.{library}.{chromosome}.log",
    benchmark:
        GATK / "apply_bqsr/{sample}.{library}.{chromosome}.bmk"
    conda:
        "../envs/gatk4.yml"
    params:
        extra=params["gatk4"]["apply_bqsr"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk ApplyBQSR \
            {params.extra} \
            --input {input.bam} \
            --reference {input.reference} \
            --bqsr-recal-file {input.table} \
            --output {output.bam} \
        2> {log} 1>&2
        """


rule gatk4_apply_bqsr_all:
    """Apply the recalibration table to all libraries and chromosomes"""
    input:
        [
            GATK / f"apply_bqsr/{sample}.{library}.{chromosome}.bam"
            for sample, library in SAMPLE_LIB
            for chromosome in CHROMOSOMES
        ],


rule gatk4_apply_bqsr_report:
    """Generate a report for all libraries and chromosomes"""
    input:
        [
            GATK / f"apply_bqsr/{sample}.{library}.{chromosome}.{report}"
            for sample, library in SAMPLE_LIB
            for chromosome in CHROMOSOMES
            for report in BAM_REPORTS
        ],


# rule gatk4_analyze_covariates:
#     shell:
#         pass


rule gatk4_haplotype_caller_one:
    """Call variants for a single library and chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        bam=GATK / "apply_bqsr/{sample}.{library}.{chromosome}.bam",
        dict_=REFERENCE / "genome.dict",
    output:
        gvcf_gz=GATK / "haplotype_caller/{sample}.{library}.{chromosome}.gvcf.gz",
    log:
        GATK / "haplotype_caller/{sample}.{library}.{chromosome}.log",
    benchmark:
        GATK / "haplotype_caller/{sample}.{library}.{chromosome}.bmk"
    conda:
        "../envs/gatk4.yml"
    params:
        extra=params["gatk4"]["haplotype_caller"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk HaplotypeCaller \
            {params.extra} \
            --reference {input.reference} \
            --input {input.bam} \
            --output {output.gvcf_gz} \
            --emit-ref-confidence GVCF \
        2> {log} 1>&2
        """


rule gatk4_haplotype_caller_all:
    """Call variants for all libraries and chromosomes"""
    input:
        [
            GATK / f"haplotype_caller/{sample}.{library}.{chromosome}.gvcf.gz"
            for sample, library in SAMPLE_LIB
            for chromosome in CHROMOSOMES
        ],


# I think here comes the GenomicsDBImport step


rule gatk4_combine_gvcfs_one:
    """Combine gVCFs to get a chromosome"""
    input:
        vcf_gzs=get_files_to_genotype,
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf_gz=GATK / "joint_variants/{chromosome}.vcf.gz",
    log:
        GATK / "joint_variants/{chromosome}.log",
    benchmark:
        GATK / "joint_variants/{chromosome}.bmk"
    conda:
        "../envs/gatk4.yml"
    params:
        variant_line=compose_v_line,
        extra=params["gatk4"]["combine_gvcfs"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk CombineGVCFs \
            {params.extra} \
            {params.variant_line} \
            --reference {input.reference} \
            --output {output.vcf_gz} \
        2> {log} 1>&2
        """


rule gatk4_combine_gvcfs_all:
    """Get all chromosomal gVCFs"""
    input:
        [GATK / f"joint_variants/{chromosome}.vcf.gz" for chromosome in CHROMOSOMES],


rule gatk4_genotype_gvcfs_one:
    """Genotype a single chromosome"""
    input:
        vcf_gz=GATK / "joint_variants/{chromosome}.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf_gz=GATK / "genotyped_variants/{chromosome}.vcf.gz",
    log:
        GATK / "genotyped_variants/{chromosome}.log",
    benchmark:
        GATK / "genotyped_variants/{chromosome}.bmk"
    conda:
        "../envs/gatk4.yml"
    params:
        extra=params["gatk4"]["genotype_gvcfs"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk GenotypeGVCFs \
            {params.extra} \
            --variant {input.vcf_gz} \
            --reference {input.reference} \
            --output {output.vcf_gz} \
        2> {log} 1>&2
        """


rule gatk4_genotype_gvcfs_all:
    """Genotype all chromosomes"""
    input:
        [GATK / f"genotyped_variants/{chromosome}.vcf.gz" for chromosome in CHROMOSOMES],


# According to GATK Best Practices, we should be using VariantRecalibrator and
# Apply VQSR, but since they need datasets with known variants (human and
# genome-wide), we skip these steps.
# According to GATK, https://gatk.broadinstitute.org/hc/en-us/articles/360035532412-Can-t-use-VQSR-on-non-model-organism-or-small-dataset
# We have to :
#   - hard filter variants instead, and/or
#   - use CNNs to learn from the data
# So far, CNNs are only fit for single-sample SNPs, so we are voiding it


rule gatk4_calculate_genotype_posteriors_one:
    """Calculate genotype posteriors for a single chromosome"""
    input:
        vcf=GATK / "genotyped_variants/{chromosome}.vcf.gz",
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
    output:
        vcf=GATK / "variants_posteriors/{chromosome}.vcf.gz",
    log:
        GATK / "variants_posteriors/{chromosome}.log",
    benchmark:
        GATK / "variants_posteriors/{chromosome}.bmk"
    conda:
        "../envs/gatk4.yml"
    params:
        extra=params["gatk4"]["calculate_genotype_posteriors"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk CalculateGenotypePosteriors \
            {params.extra} \
            --output {output.vcf} \
            --variant {input.vcf} \
            --reference {input.reference} \
        2> {log} 1>&2
        """


rule gatk4_calculate_genotype_posteriors_all:
    """Calculate genotype posteriors for all chromosomes"""
    input:
        [
            GATK / f"variants_posteriors/{chromosome}.vcf.gz"
            for chromosome in CHROMOSOMES
        ],


rule gatk4_variant_filtration_one:
    """Filter variants for a single chromosome"""
    input:
        reference=REFERENCE / "genome.fa.gz",
        dict_=REFERENCE / "genome.dict",
        vcf=GATK / "variants_posteriors/{chromosome}.vcf.gz",
    output:
        vcf=GATK / "variants_filtered/{chromosome}.vcf.gz",
    log:
        GATK / "variants_filtered/{chromosome}.log",
    benchmark:
        GATK / "variants_filtered/{chromosome}.bmk"
    conda:
        "../envs/gatk4.yml"
    params:
        filter_name=params["gatk4"]["variant_filtration"]["filter_name"],
        filter_expression=params["gatk4"]["variant_filtration"]["filter_expression"],
        extra=params["gatk4"]["variant_filtration"]["extra"],
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        """
        gatk VariantFiltration \
            {params.extra} \
            --reference {input.reference} \
            --variant {input.vcf} \
            --output {output.vcf} \
            --filter-expression '{params.filter_expression}' \
            --filter-name '{params.filter_name}' \
        2> {log} 1>&2
        """


rule gatk4_variant_filtration_all:
    """Filter variants for all chromosomes"""
    input:
        [GATK / f"variants_filtered/{chromosome}.vcf.gz" for chromosome in CHROMOSOMES],


# Now should come funcotator, but it is human centric
# We use snpeff instead
# rule gatk4_funcotator


# rule gatk_cnn_score_variants:
#     input:
#         bams = [
#             GATK / f"{sample}.{library}.bqsr.bam"
#             for sample, library in SAMPLE_LIB
#         ],
#         bais = [
#             GATK / f"{sample}.{library}.bqsr.bam.bai"
#             for sample, library in SAMPLE_LIB
#         ],
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


rule gatk4_variant_filtration_merge:
    """Merge all VCF chromosomes"""
    input:
        expand(
            GATK / "variants_filtered/{chromosome}.vcf.gz",
            chromosome=CHROMOSOMES,
        ),
    output:
        GATK / "variants_filtered.vcf.gz",
    log:
        GATK / "variants_filtered.log",
    conda:
        "../envs/bcftools.yml"
    threads: 24
    shell:
        """
        bcftools concat \
            --naive \
            --output {output} \
            --output-type z9 \
            --threads {threads} \
            {input} \
        2> {log} 1>&2
        """


rule gatk4_all:
    input:
        GATK / "variants_filtered.vcf.gz",


rule gatk4_report:
    input:
        rules.gatk4_base_recalibrator_all.input,
        rules.gatk4_apply_bqsr_report.input,


rule gatk4:
    """Run all GATK4 steps"""
    input:
        # rules.gatk4_base_recalibrator_all.input,
        # rules.gatk4_apply_bqsr_all.input,
        # rules.gatk4_haplotype_caller_all.input,
        # rules.gatk4_combine_gvcfs.output,
        # rules.gatk4_genotype_gvcfs.output,
        # GATK / "genotyped_variants.vcf.gz",
        # GATK / "variants_posteriors.vcf.gz",
        GATK / "variants_filtered.vcf.gz",
