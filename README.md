# Snakemake workflow: `Bioinfo_Macro_Host_Genomics`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/3d-omics/Bioinfo_Macro_Host_Genomics/workflows/Tests/badge.svg?branch=main)](https://github.com/3d-omics/Bioinfo_Macro_Host_Genomics/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for  Short Variant Discovery in Host Genomes


## Usage

- Test that it works:
  - Make sure you have installed snakemake, samtools and bcftools. Either
    - install them with conda/mamba :`conda install -c bioconda samtools bcftools`).
    - or create an environment (`conda create -n 3dohg -c bioconda snakemake samtools bcftools`), and activate it (`conda activate 3dohg`)
  - Generate mock data with `bash workflow/scripts/generate_mock_data.sh`
  - Run the pipeline: `snakemake --use-conda --jobs 8 all`. It will download all the necesary software through conda. It should take less than 5 minutes.

- Run it with your own data:
  - Edit `config/samples.tsv` and add your samples and where are they located.
  - Edit `config/features.tsv` with information regarding the reference you are using.
  - Run the pipeline: `snakemake --use-conda --jobs 8 all`.
  - (slurm users): `./run_slurm`

## Features

- FASTQ processing with [`fastp`](https://github.com/OpenGene/fastp)
- Mapping with [`bowtie2`](https://github.com/BenLangmead/bowtie2)
- SAM/BAM/CRAM processing with [`samtools`](https://github.com/samtools/samtools) and [`picard`](https://github.com/broadinstitute/picard)
- Sample swap detection with [`gtcheck`](https://github.com/samtools/bcftools)
- SNP calling with [`GATK4`](https://github.com/broadinstitute/gatk)
- SNP annotation with [`SNPEff`](https://github.com/pcingola/SnpEff)

## DAG

![host_genomics_pipeline](./rulegraph.svg?raw=true)

## References

TBA
