rule reads_link_one:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=READS / "{sample}.{library}_1.fq.gz",
        reverse_=READS / "{sample}.{library}_2.fq.gz",
    log:
        READS / "{sample}.{library}.log",
    conda:
        "../envs/empty.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_}
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_}
        """


rule reads_link_all:
    """Link all reads in the samples.tsv"""
    input:
        [
            READS / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
        ],


rule reads_fastqc_all:
    """Run fastqc on all raw reads"""
    input:
        [
            READS / f"{sample}.{library}_{end}_fastqc.{extension}"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads_link_all.input,
        rules.reads_fastqc_all.input,


localrules:
    reads_link_one,
