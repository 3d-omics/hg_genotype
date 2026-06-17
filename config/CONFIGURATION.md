# Configuration guide

Three files must be filled before running the pipeline:

- [`config/samples.tsv`](#samplestsv) — one row per library
- [`config/features.yml`](#featuresyml) — reference genome and chromosome metadata
- [`config/params.yml`](#paramsyml) — tool parameters

---

## samples.tsv

A tab-separated file. Each row is one sequencing library. A sample with multiple
libraries (e.g. re-sequenced) gets one row per library; they are merged at the
mark-duplicates step.

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Unique sample identifier |
| `library_id` | yes | Library identifier (unique within a sample) |
| `forward_filename` | yes | Path to R1 FASTQ (`.fq.gz`) |
| `reverse_filename` | yes | Path to R2 FASTQ (`.fq.gz`) |
| `sex` | yes | Sex label — must match a key in `sex_ploidy` (see features.yml), or any string if `sex_ploidy` is empty |
| `ploidy` | yes | Per-individual ploidy for autosomes (usually `2`) |
| `pool_size` | yes | Number of pooled individuals (usually `1`; set >1 for pool-seq) |

Lines starting with `#` are ignored.

### Checklist

- [ ] Every sample has at least one library row
- [ ] `forward_filename` and `reverse_filename` paths exist and are readable
- [ ] `sex` values match the keys defined in `features.yml` → `sex_ploidy` (or `sex_ploidy` is omitted)
- [ ] `ploidy` is set to the individual's autosomal ploidy (2 for diploid, 1 for haploid, 4 for tetraploid, …)
- [ ] `pool_size` is `1` for individual samples, or the pool size for pool-seq experiments

### Examples

**Standard diploid samples, ZW sex system (e.g. birds):**

```tsv
sample_id	library_id	forward_filename	reverse_filename	sex	ploidy	pool_size
hen1	lib1	resources/reads/hen1_1.fq.gz	resources/reads/hen1_2.fq.gz	female	2	1
rooster1	lib1	resources/reads/rooster1_1.fq.gz	resources/reads/rooster1_2.fq.gz	male	2	1
```

**Sample sequenced across two lanes (multiple libraries):**

```tsv
sample_id	library_id	forward_filename	reverse_filename	sex	ploidy	pool_size
sample1	lib1	resources/reads/sample1_lane1_1.fq.gz	resources/reads/sample1_lane1_2.fq.gz	female	2	1
sample1	lib2	resources/reads/sample1_lane2_1.fq.gz	resources/reads/sample1_lane2_2.fq.gz	female	2	1
```

**Pool-seq (10 diploid individuals per pool):**

```tsv
sample_id	library_id	forward_filename	reverse_filename	sex	ploidy	pool_size
pool_A	lib1	resources/reads/poolA_1.fq.gz	resources/reads/poolA_2.fq.gz	female	2	10
```

**Haploid organism with no sex chromosomes:**

```tsv
sample_id	library_id	forward_filename	reverse_filename	sex	ploidy	pool_size
strain1	lib1	resources/reads/strain1_1.fq.gz	resources/reads/strain1_2.fq.gz	unknown	1	1
```

---

## features.yml

Describes the reference genome and how chromosomes should be treated.

### Fields

| Field | Required | Description |
|---|---|---|
| `name` | yes | Short assembly name; used as a filename prefix throughout `results/` |
| `species` | yes | Lowercase species name for the VEP cache download (e.g. `homo_sapiens`) |
| `release` | yes | Ensembl release number for the VEP cache (e.g. `110`) |
| `dna` | yes | Path to the reference FASTA (`.fa.gz`) |
| `gtf` | yes | Path to the annotation GTF (`.gtf.gz`). Must be GTF, not GFF3 |
| `known_vcf` | yes | Path to a VCF of known variants for BQSR (`.vcf.gz`). Can be bootstrapped with `results/align/bcftools/bcftools.vcf.gz` |
| `regions` | yes | Path to a BED4 file listing all regions to genotype (one region per row, 4th column = region name) |
| `organelles` | yes | List of organellar chromosome names (mitochondria, chloroplast). These are called at `pool_size` ploidy regardless of individual ploidy |
| `sex_ploidy` | no | Per-sex ploidy table for sex chromosomes. Any chromosome absent from a sex's entry is mock-called (ploidy 0). Omit entirely for organisms with no sex chromosomes |

### Checklist

- [ ] `name` contains no spaces or special characters (it becomes a filename prefix)
- [ ] `species` matches the Ensembl species slug exactly (check at ensembl.org)
- [ ] `release` matches an available Ensembl release for that species
- [ ] `dna`, `gtf`, `known_vcf`, `regions` paths exist and are readable
- [ ] `gtf` is a GTF file, not GFF3
- [ ] `organelles` lists all non-nuclear chromosomes present in the BED4
- [ ] All sex labels used in `samples.tsv` appear as keys in `sex_ploidy` (if `sex_ploidy` is set)
- [ ] Every sex chromosome present in the BED4 appears in at least one sex's `sex_ploidy` entry

### Examples

**Birds (ZW system):**

```yaml
name: GRCg6a
species: gallus_gallus
release: 106
dna: resources/reference/GRCg6a.fa.gz
gtf: resources/reference/GRCg6a.gtf.gz
known_vcf: resources/reference/GRCg6a.known.vcf.gz
regions: resources/reference/GRCg6a.bed4
organelles: [MT]
sex_ploidy:
  male:
    Z: 2
  female:
    Z: 1
    W: 1
```

**Mammals (XY system):**

```yaml
name: GRCh38
species: homo_sapiens
release: 110
dna: resources/reference/GRCh38.fa.gz
gtf: resources/reference/GRCh38.gtf.gz
known_vcf: resources/reference/GRCh38.known.vcf.gz
regions: resources/reference/GRCh38.bed4
organelles: [MT]
sex_ploidy:
  male:
    X: 1
    Y: 1
  female:
    X: 2
```

**X0 system (no Y chromosome in reference):**

```yaml
sex_ploidy:
  male:
    X: 1
  female:
    X: 2
```

**Haploid organism, no sex chromosomes:**

```yaml
name: MyAssembly_v1
species: my_species
release: 100
dna: resources/reference/assembly.fa.gz
gtf: resources/reference/assembly.gtf.gz
known_vcf: resources/reference/assembly.known.vcf.gz
regions: resources/reference/assembly.bed4
organelles: [MT]
```

**Plant (with chloroplast):**

```yaml
organelles: [MT, Pt]
```

---

## params.yml

Controls filtering thresholds and tool extra arguments. Defaults are suitable for
most projects; only change them if you have a specific reason.

### Fields

| Field | Description |
|---|---|
| `align.bwamem2.extra` | Extra flags passed to `bwa-mem2 mem` |
| `variants.filter.SNP` | GATK hard-filter expression for SNPs |
| `variants.filter.INDEL` | GATK hard-filter expression for INDELs |

### Checklist

- [ ] GATK filter expressions use valid JEXL syntax
- [ ] `bwamem2.extra` does not duplicate flags already set by the wrapper (e.g. `-R` for read groups is added automatically)

### Example

```yaml
align:
  bwamem2:
    extra: ""

variants:
  filter:
    SNP: "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
    INDEL: "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
```
