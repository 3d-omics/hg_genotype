#!/usr/bin/env bash
set -euo pipefail

mkdir --parents resources/reads/ resources/reference/

# Avoid using the telomeres
genome_region="1:490000-500000"
reads_per_sample=3000

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
gzip \
    --decompress \
    --stdout \
    Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.1.fa.gz \
| bgzip --threads 8 \
> ggal.fa.gz

gzip \
    --decompress \
    --stdout \
    gallus_gallus.vcf.gz \
| bgzip --threads 8 \
> ggal.vcf.gz

# Index
samtools faidx ggal.fa.gz
bcftools index ggal.vcf.gz

# Slice the 1st Mbp from chromosome 1
samtools faidx ggal.fa.gz "$genome_region" \
| sed "s/$genome_region/1/g" \
| bgzip --threads 8 \
> ggal.mock.fa.gz

bcftools view ggal.vcf.gz "$genome_region" \
| bgzip --threads 8 \
> ggal.mock.vcf.gz


popd


# Simulate reads with wgsim. Different seed means different individuals
for individual in sample{1..2}; do

    wgsim \
        -N "$reads_per_sample" \
        -1 100 \
        -2 100 \
        -S "$individual" \
        resources/reference/ggal.mock.fa.gz \
        >(pigz -9 > resources/reads/"${individual}"_1.fq.gz) \
        >(pigz -9 > resources/reads/"${individual}"_2.fq.gz) \
    | pigz -9 > resources/reads/"${individual}".tsv.gz

done
