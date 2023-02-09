#!/usr/bin/env bash
set -euo pipefail

mkdir --parents resources/reads/ resources/reference/

pushd resources/reference/
# Download reference
wget \
    --continue \
    https://ftp.ensembl.org/pub/release-108/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.1.fa.gz


# known variants
wget \
    --continue \
    https://ftp.ensembl.org/pub/release-108/variation/vcf/gallus_gallus/gallus_gallus.vcf.gz


# Recompress with bgzip
gzip --decompress --stdout Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.1.fa.gz | bgzip --threads 8 > ggal.fa.gz
gzip --decompress --stdout gallus_gallus.vcf.gz | bgzip --threads 8 > ggal.vcf.gz

# Index
samtools faidx ggal.fa.gz
bcftools index ggal.vcf.gz

# Slice the 1st Mbp from chromosome 1
samtools faidx ggal.fa.gz 1:1-1000000 | sed "s/1:1-1000000/1/g" | bgzip --threads 8 > ggal.mock.fa.gz
bcftools view ggal.vcf.gz 1:1-1000000 | bgzip --threads 8 > ggal.mock.vcf.gz


popd


# Simulate reads with wgsim. Different seed means different individuals

for individual in sample{1..2}; do
    wgsim \
        -N 100000 \
        -1 150 \
        -2 150 \
        -S "$individual" \
        resources/reference/ggal.mock.fa.gz \
        >(pigz -9 > resources/reads/"${individual}"_1.fq.gz) \
        >(pigz -9 > resources/reads/"${individual}"_2.fq.gz) \
    | pigz -9 > resources/reads/"${individual}".tsv.gz
done
