#!/usr/bin/env bash
set -euo pipefail

mkdir --parents resources/reads/ resources/reference/

# Avoid using the telomeres
reads_per_sample=30000

pushd resources/reference/
# Download reference
wget \
    --continue \
    https://ftp.ensembl.org/pub/release-111/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz


# known variants
wget \
    --continue \
    https://ftp.ensembl.org/pub/release-111/variation/vcf/gallus_gallus/gallus_gallus.vcf.gz


# Recompress with bgzip
gzip \
    --decompress \
    --stdout \
    Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz \
| bgzip \
    --threads 8 \
> ggal.fa.gz

gzip \
    --decompress \
    --stdout \
    gallus_gallus.vcf.gz \
| bgzip \
    --threads 8 \
> ggal.vcf.gz


# Index
samtools faidx ggal.fa.gz
bcftools index ggal.vcf.gz


# Slice the dataset
samtools faidx \
    ggal.fa.gz \
    1:490000-500000 \
    2:490000-500000 \
    Z:490000-500000 \
    W:490000-500000 \
    MT \
| seqtk seq \
| cut \
    -f 1 \
    -d : \
| bgzip \
    --threads 8 \
> ggal.mock.fa.gz

bcftools view \
    ggal.vcf.gz \
    1:490000-500000 \
    2:490000-500000 \
    Z:490000-500000 \
    W:490000-500000 \
    MT \
| awk '{
    if ($0 ~ /^#|^MT/) {print $0}
    else{$2 = $2 - 490000; print $0}
    }' \
    OFS="\t" \
| bgzip \
    --threads 8 \
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
