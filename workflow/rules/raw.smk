rule raw_link_pe:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=READS + "/{sample}.{library}_1.fq.gz",
        reverse_=READS + "/{sample}.{library}_2.fq.gz",
    log:
        READS + "{sample}.{library}.log",
    conda:
        CONDA + "empty.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_}
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_}
        """


rule raw_link:
    input:
        [
            f"{READS}/{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
        ],


rule raw_fastqc:
    input:
        [
            f"{READS}/{sample}.{library}_{end}_fastqc.html"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
        ],


# rule raw_extract_genome:
#     """Extract the fasta.gz on config.yaml into genome.fa"""
#     input:
#         fa_gz = features["reference"]["dna"]
#     output:
#         fa = READS + "genome.fa"
#     log:
#         RAW + "genome.log"
#     benchmark:
#         RAW + "genome.json"
#     shell:
#         "gzip --decompress --stdout {input.fa_gz} > {output.fa} 2> {log}"
