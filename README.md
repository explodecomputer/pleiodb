# pleiodb

Compact storage and query engine for GWAS summary statistics across many variants and traits.
Designed for V ≈ 100,000 variants × T ≈ 20,000 traits; stores z-scores, effective sample size,
and effect allele frequency in chunk-compressed binary format (~5–6 GB on disk at full scale).

## Installation

Requires [mamba](https://mamba.readthedocs.io) (or conda). `bcftools` and `cyvcf2` must come from bioconda.

```bash
mamba env create -f environment.yml
mamba activate pleiodb
```

## Input files

**Variants TSV** — two columns, no header: `ALID` and `EAF`.
ALID format: `CHROM:POS_A1_A2` where A1 ≤ A2 alphabetically (e.g. `10:101558746_G_T`).

```
10:101558746_G_T	0.56079
10:1206798_C_T	0.071479
19:44908822_C_T	0.12
```

**Traits TSV** — four columns, no header: `trait_id`, `trait_name`, `vcf_path`, `build` (optional).
VCF files must follow the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification).

```
ieu-b-2	Body mass index	/data/gwas/bmi.vcf.gz	hg19
ukb-b-1234	LDL cholesterol	/data/gwas/ldl.vcf.gz	hg38
finn-r-T2D		/data/gwas/finngen_T2D.vcf.gz
```

## Quick start

The repo includes test data (501 variants, 5 traits) for a self-contained walkthrough.

```bash
# Build a database
pleiodb build my.pleiodb \
  --variants tests/test_data/variants_hg19.tsv \
  --traits   tests/test_data/traits.tsv

# Inspect the database
pleiodb info my.pleiodb

# All traits for a single variant
pleiodb query my.pleiodb --variant 10:1206798_C_T

# All associations passing p ≤ 1×10⁻⁵
pleiodb query my.pleiodb --pval 1e-5

# Single trait, genome-wide
pleiodb query my.pleiodb --trait ieu-a-7
```

## CLI reference

```
pleiodb build   OUTPUT  --variants TSV  --traits TSV  [--workers N]  [--variants-build BUILD]
pleiodb query   DB      --variant ID | --trait ID | --region CHR:START-END | --pval FLOAT
                        [--variants-file FILE]  [--traits-file FILE]
                        [--beta-se]  [--format tsv|json]  [-o FILE]
pleiodb lambda  DB      [--null-thresh FLOAT]  [--workers N]
pleiodb info    DB
```

Pass `--variants-build hg38` (or `hg19`) when your variant list and VCF files are on different
genome builds — pleiodb will lift coordinates automatically using `pyliftover`.

## Running the tests

```bash
pytest tests/test_integration.py -v
```
