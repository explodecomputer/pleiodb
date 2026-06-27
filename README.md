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

**Traits TSV** — tab-separated with header row.
Required columns: `trait_id`, `trait_name`, `N`, `vcf_path`.
Optional columns: `K` (case fraction, required for binary traits), `build`.
VCF files must follow the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification).

```
trait_id    trait_name          N       K     vcf_path                    build
ieu-b-2     Body mass index     461460        /data/gwas/bmi.vcf.gz       hg19
ukb-b-1234  LDL cholesterol     94595         /data/gwas/ldl.vcf.gz       hg19
finn-r-T2D  Type 2 diabetes     200000  0.15  /data/gwas/finngen_T2D.vcf.gz
```

`K` is the case fraction (`N_cases / N`); omit or leave blank for continuous traits.
Trait IDs must not contain commas (commas are used as the delimiter for `--traits`).

## Quick start

The repo includes test data (501 variants, 5 traits) for a self-contained walkthrough.

```bash
# Build a database
pleiodb build my.pleiodb \
  --variants tests/test_data/variants_hg19.tsv \
  --traits   tests/test_data/traits.tsv

# Inspect the database (warns if rho matrix is absent)
pleiodb info my.pleiodb

# Compute the rho matrix (sample-overlap-weighted phenotypic correlation)
pleiodb rho my.pleiodb

# All traits for a single variant
pleiodb query my.pleiodb --variant 10:1206798_C_T

# All associations passing p ≤ 1×10⁻⁵
pleiodb query my.pleiodb --pval 1e-5

# Single trait, genome-wide
pleiodb query my.pleiodb --trait ieu-a-7

# All associations passing p ≤ 1×10⁻⁵ for single trait and variant
pleiodb query my.pleiodb --pval 1e-5 --trait ukb-b-10787 --variant 2:25121853_A_G

# Intersect for multiple traits and variants
echo -e "6:36454223_A_T\n1:89472196_A_G" > v.txt
echo -e "ukb-b-10787\nieu-a-7" > t.txt
pleiodb query my.pleiodb --traits-file t.txt --variants-file v.txt

# Query rho for a specific pair of traits
pleiodb rho my.pleiodb --traits ieu-a-7,ukb-b-10787

# Query rho for a set of traits, output as a square matrix
pleiodb rho my.pleiodb --traits-file t.txt --matrix

# Query rho as JSON
pleiodb rho my.pleiodb --traits ieu-a-7,ukb-b-10787 --format json
```

## The rho matrix

`pleiodb rho` computes and stores a symmetric T×T matrix where each cell
`rho[j, k]` estimates `ρ_pheno × N_overlap[j,k] / √(N_j × N_k)` — the
phenotypic correlation between traits j and k weighted by their shared sample
fraction. This is used downstream for z-score correction, multivariate models,
pleiotropy tests (GSMR, MR-RAPS), and conditional F-statistics.

**Estimator**: conditional maximum likelihood (Forde et al., `mr.simss::est_lambda`),
which corrects the truncation bias of a naive Pearson correlation on null
z-scores. Null variants are those with |z_j| < 1.0 AND |z_k| < 1.0 (default
threshold; adjust with `--null-thresh`). Pairs with fewer than 500 shared null
variants return NaN (override with `--min-nulls`).

```bash
# Compute with non-default settings
pleiodb rho my.pleiodb --null-thresh 0.5 --min-nulls 200 --workers 16

# Query: pairwise list (default)
pleiodb rho my.pleiodb --traits ieu-a-7,ukb-b-10787,ukb-b-19953
# trait_id_1    trait_id_2    rho
# ieu-a-7       ukb-b-10787   -0.3967
# ieu-a-7       ukb-b-19953   -0.0995
# ukb-b-10787   ukb-b-19953    0.1125

# Query: square matrix
pleiodb rho my.pleiodb --traits ieu-a-7,ukb-b-10787,ukb-b-19953 --matrix
```

`pleiodb info` reports `rho_present` and emits a warning when the rho matrix
has not yet been computed.

## LD-based z-score imputation

When `--ld-dir` is supplied at build time, pleiodb runs a post-build imputation pass that
fills missing z-scores using elastic-net regression on LD eigenvectors.

**How it works** (dense VCF mode):

1. After all traits are ingested, for each LD block pleiodb identifies traits that have
   missing z-scores at positions covered by that block.
2. For those traits it reads z-scores for **all** variants in the LD block directly from
   the source GWAS-VCF files — not just the ~95k stored variants, but the full set of
   reference panel variants (~5–10k per block).
3. Precomputed block eigenvectors (stored as `.ldeig.rds` alongside the LD panel) are
   loaded and cached as `.ldeig.npz` on first use (~62 ms cached vs ~3–5 s cold).
4. An elastic-net model is trained on the dense z-score vector using the eigenvectors
   as features; predictions are extracted only for the missing stored positions.
5. A polynomial rescaling step corrects bias and variance of the imputed values.

**Why dense beats sparse**: the stored variant set covers only 20–50 of the ~5–10k
reference panel SNPs per block (0.4–1%). Training on the full dense set gives the
elastic net far more signal to learn the local LD pattern.

**Benchmark (chr1, 5 traits, 32 workers):** 0.71 s/block with dense VCF mode.

```bash
# Build with LD imputation (requires LD reference panel and hg38 variant list)
pleiodb build my.pleiodb \
  --variants variants_hg38.tsv \
  --traits   traits.tsv \
  --ld-dir   /data/ld_reference_panel_hg38/EUR \
  --workers  32 \
  --ld-vcf-threads 8
```

The traits TSV must include a `vcf_path` column (and optionally a `build` column for
hg19 VCFs) so that pleiodb can read back the source files during imputation:

```
trait_id    trait_name          N       K     vcf_path                    build
ieu-b-2     Body mass index     461460        /data/gwas/bmi.vcf.gz       hg38
ukb-b-1234  LDL cholesterol     94595         /data/gwas/ldl.vcf.gz       hg19
```

Imputed positions are recorded in `imputed.coo.zst` inside the database directory
(COO-format sparse array of (variant_idx, trait_idx) pairs). `pleiodb info` reports
the total imputed cell count.

## CLI reference

```
pleiodb build   OUTPUT  --variants TSV  --traits TSV  [--workers N]  [--variants-build BUILD]
                        [--ld-dir DIR]  [--ld-vcf-threads N]  [--ld-thresh F]  [--ld-min-cor F]
pleiodb rho     DB      [--null-thresh FLOAT]  [--min-nulls INT]  [--workers N]
                        [--traits t1,t2,...]   [--traits-file FILE]
                        [--matrix]  [--format tsv|json]  [-o FILE]
pleiodb query   DB      --variant ID | --trait ID | --region CHR:START-END | --pval FLOAT
                        [--variants-file FILE]  [--traits-file FILE]
                        [--format tsv|json]  [-o FILE]
pleiodb info    DB
```

Pass `--variants-build hg38` (or `hg19`) when your variant list and VCF files are on different
genome builds — pleiodb will lift coordinates automatically using `pyliftover`.

## Running the tests

```bash
pytest tests/test_integration.py -v
```
