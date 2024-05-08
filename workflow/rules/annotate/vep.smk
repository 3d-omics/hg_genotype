rule annotate__vep__:
    input:
        vcf=FILTER / "all.filtered.vcf.gz",
        fa=REFERENCE / "genome.fa.gz",
        gtf=REFERENCE / "annotation.gtf.gz",
        gtf_tbi=REFERENCE / "annotation.gtf.gz.tbi",
    output:
        vcf=VEP / "{sample}.vcf.gz",
        tmp=temp(VEP / "{sample}.tmp.vcf.gz"),
        html=VEP / "{sample}.vep.html",
    log:
        VEP / "{sample}.log",
    conda:
        "__environment__.yml"
    params:
        sample=lambda w: w.sample,
    shell:
        """
        bcftools view \
            --samples {params.sample} \
            --output-type z \
            --output-file {output.tmp} \
            --trim-alt-alleles \
            {input.vcf} \
        2> {log} 2>&1

        vep \
            --input_file {output.tmp} \
            --output_file {output.vcf} \
            --fasta {input.fa} \
            --gtf {input.gtf} \
            --compress_output bgzip \
            --stats_file {output.html} \
            --buffer_size 500 \
        2>> {log} 1>&2
        """


rule annotate__vep__reports:
    input:
        [VEP / f"{sample}.vep.html" for sample in SAMPLES],


rule annotate__vep:
    input:
        [VEP / f"{sample}.vcf.gz" for sample in SAMPLES],
