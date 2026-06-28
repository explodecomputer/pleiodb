# pleiodb

A compact binary store and query engine for GWAS summary statistics across many variants and traits simultaneously.

## Language

### Variants and alleles

**Variant**:
A biallelic locus identified by `CHROM:POS_REF_ALT` where REF and ALT are in canonical allele orientation. The identifier encodes all coordinates; no external reference genome is required.
_Avoid_: SNP (too narrow — indels are also stored), site

**Canonical allele orientation**:
The convention where REF = `min(ref, alt)` alphabetically and ALT = `max(ref, alt)`. When a source VCF has the alleles in non-canonical orientation (genome-reference REF > ALT alphabetically), the z-score is negated on ingestion to compensate.
_Avoid_: "reference-aligned", "forward strand"

**ALID (Alphabetic Variant ID)**:
The string `CHROM:POS_A1_A2` where A1 < A2 alphabetically (e.g. `10:101558746_G_T`). The canonical primary key for a variant in pleiodb. Encodes chromosome, position, and both alleles; all coordinates are parseable from the ALID alone. For alleles longer than 20 characters, the allele is compressed to `{allele[:8]}~{sha256(allele)[:4]}` (13 chars) — the tilde is the marker. Collision resistance is only required within a single genomic position, where at most a handful of long alleles exist, making 4 hex chars (16 bits) sufficient.
_Avoid_: rsID (rsIDs are an alternate lookup key, not the ALID), "variant ID" (use ALID to be precise)

**Positional matching**:
The strategy for aligning variant IDs to VCF records: match by CHROM:POS first (fast index lookup), then confirm by REF/ALT alleles. rsIDs in the VCF `ID` field are not used for primary matching.
_Avoid_: ID-based matching

### Traits and data

**Trait**:
A GWAS phenotype represented as a column in the V×T matrix. Identified by a `trait_id` string (e.g. `ieu-a-7`). Trait IDs must not contain commas — commas are used as the delimiter when supplying a list of trait IDs on the command line.

**GWAS-VCF**:
A VCF file following the MRC IEU GWAS-VCF specification. Source of z-scores and sample sizes per variant per trait.

**Z-score**:
The primary stored statistic: effect size (ES) divided by standard error (SE), or the pre-computed `EZ` field when present.

**Trait type**:
Derived from the traits input file: if `K` (case fraction) is present, the trait is binary (logistic regression, betas on log-OR scale); if absent, it is continuous. Never inferred from VCF fields.

**Neff (effective sample size)**:
Per-variant, per-trait value stored as a log2-encoded uint16. Derived at ingest as `var_y[t] / (SE[v,t]² × 2·EAF[v]·(1−EAF[v]))`. For a well-specified study this equals the actual sample size (continuous) or `4·N·K·(1−K)` (binary), modulated by per-variant SE variation. Not read from the VCF `SS` field.
_Avoid_: treating Neff as raw sample size — it is the effective N on the normalised (var(y)=1) scale, consistent with `SE_norm = 1/sqrt(Neff × 2p(1−p))`.

**Normalised beta**:
Beta (and SE) reconstructed on the var(y)=1 scale: `SE_norm = 1/sqrt(Neff × 2p(1−p))`, `beta_norm = z × SE_norm`. The default query output alongside Z and p-value. For binary traits, this is on an approximately normalised log-OR scale; var(y) is not explicitly set to 1 but emerges from the normalisation.
_Avoid_: "standardised beta" (ambiguous — sometimes means divided by phenotypic SD, sometimes by genotypic SD)

**Study-scale beta**:
Beta (and SE) on the original GWAS units: `beta_study = sqrt(var_y) × beta_norm`. Returned only when explicitly requested via a CLI flag. Not part of the default output.
_Avoid_: "original beta" (original to what?), "unstandardised beta"

**Default query output**:
Every query returns Z, normalised beta, normalised SE, and p-value (two-sided, derived as `2·Φ(−|Z|)`). Study-scale beta is an opt-in via a CLI flag and requires `var_y` to be present in the database.

**var_y (per-study phenotypic variance)**:
Estimated at ingest as `median_v(SE[v,t]² × 2·EAF[v]·(1−EAF[v]) × Neff_study[t])` across variants with EAF in (0.01, 0.99) and valid SE, where `Neff_study[t]` = `N[t]` for continuous traits and `4·N[t]·K[t]·(1−K[t])` for binary traits (N and K from traits input file). The `p(1−p)` term is symmetric so A1 and A2 frequency give the same numerical result — but EAF is used explicitly because pleiodb always tracks effects with respect to A2 (the effect/ALT allele). This symmetry applies only to the variance calculation; the sign of z-scores and betas is not symmetric and is always defined with respect to A2 allele count. Stored as a T-length float32 array. Required to reconstruct study-scale beta. Computed identically for continuous and binary traits — for binary traits (logistic regression) the estimate approximates `K(1−K)` where K is prevalence, and study-scale beta is on the log-OR scale. Normalised beta treats var(y)=1 for all trait types.

**Ingestion requirements**: `ES` and `SE` FORMAT fields must be present in a GWAS-VCF; a trait whose VCF lacks either is rejected at build time. `ES`/`SE` are always extracted at ingest regardless of whether `EZ` is also present — `EZ` is used for the stored z-score (preferred for precision), but `ES`/`SE` are used to derive Neff and var_y. `N` must be present in the traits input file; for binary traits, `K` (case fraction) must also be present. The VCF `SS` field is ignored.

**EAF (effect allele frequency)**:
Allele frequency for the effect allele (A2 — the alphabetically second allele). The primary per-variant frequency input; stored as float16. Used for beta/SE reconstruction. When ingesting a VCF: if the VCF alleles are already in canonical orientation (A1=REF < ALT=A2), EAF = VCF `AF`; if they are flipped (REF > ALT alphabetically), EAF = 1 − VCF `AF`.
_Avoid_: RAF (retired — "reference" is ambiguous between genome-ref and canonical-ref)

### Storage

**Database** (`*.pleiodb`):
The output directory containing all binary matrices, metadata, and index files for one collection of variants × traits. During an active build a `build_checkpoint.json` file is also present; it is deleted on successful completion and can be used to resume an interrupted build via `pleiodb build --resume`.

**Chunk**:
The unit of compressed storage: a rectangular submatrix of fixed size (default 512×512 cells) stored as a single Zstandard-compressed blob.

**Significance mask**:
A Zstandard-compressed COO-format sparse array of (v_idx, t_idx) pairs for all associations passing a p-value threshold. Enables O(hits) p-value queries.

**rho (sample-overlap-weighted phenotypic correlation)**:
A symmetric T×T matrix. Each cell `rho[j, k]` estimates the quantity `ρ_pheno × N_overlap[j,k] / √(N_j × N_k)` — the product of the phenotypic correlation between traits j and k and the fraction of shared samples (normalised by study sizes). This is a single combined quantity; phenotypic correlation and sample overlap cannot be separated from it without external information. Estimated from z-scores of null variants (|z_j| < z_threshold AND |z_k| < z_threshold, default z_threshold = 1.0) in the stored z-score matrix using a conditional maximum-likelihood estimator that corrects for the selection bias introduced by thresholding. Stored as `rho.bin/.cidx` (float16, symmetric). Computed as a separate post-build step via `pleiodb rho`. Used downstream for z-score correction, multivariate models, pleiotropy tests, and conditional F-statistics.
_Avoid_: "lambda" (collides with the per-trait genomic inflation factor λ = median(χ²)/0.456); "cross-trait intercept" (implies LDSC methodology specifically); "sample overlap matrix" (ignores the phenotypic correlation component); "Pearson correlation of null z-scores" (the stored value is the CML estimate, not the raw Pearson r).

**rho query interface**:
`pleiodb rho <db> --traits t1,t2,t3` or `pleiodb rho <db> --traits-file path` selects a subset of traits and returns pairwise rho values. Default output is a **pairwise list** (columns: `trait_id_1`, `trait_id_2`, `rho`) with one row per unique unordered pair. The `--matrix` flag pivots to a square matrix with trait IDs as both row labels and column headers. The `--traits`/`--traits-file` flags put the command into query mode; omitting both triggers computation of the full rho matrix (the post-build step).

## Input file formats

**Variant input file**: Two columns, tab-separated, no header: `ALID | EAF`. Provides both the variant list and per-variant effect allele frequencies in a single file. ALIDs with alleles longer than 20 characters must use the compressed form `{allele[:8]}~{sha256(allele)[:4]}`.

**Trait input file**: Tab-separated, with header. Columns: `trait_id | trait_name | N | K | vcf_path | build(optional)`. `N` is total sample size (required). `K` is case fraction `Ncase/N` (required for binary traits; if absent the trait is treated as continuous). `build` is optional genome build tag.

## Database per-trait metadata file (`traits.tsv` inside `.pleiodb`)

A human-readable tab-separated file with one row per trait. Row order defines the trait index used by the V×T binary matrices. Replaces the binary `traits.npy` and `neff_base.f32`. Columns:

| column | description |
|--------|-------------|
| `trait_id` | trait identifier |
| `trait_name` | human-readable label |
| `N` | total sample size (from input) |
| `K` | case fraction (empty for continuous) |
| `neff_study` | `N` (continuous) or `4·N·K·(1−K)` (binary) |
| `var_y` | estimated phenotypic variance (median across common-MAF variants) |
| `n_variants` | observed (pre-imputation) variants with a valid z-score for this trait |
| `n_variants_imputed` | cells filled by LD imputation (0 when `--ld-dir` is not supplied) |
| `n_variants_var_y` | variants used in the var_y median estimate (EAF ∈ (0.01, 0.99), valid SE) |
| `n_sig_{threshold}` | number of significant variants at each p-value threshold (one column per threshold, e.g. `n_sig_5e-8`, `n_sig_1e-5`) |

## Test data (`tests/test_data/`)

- **`variants_hg19.tsv`** — 500 variants in ALID+EAF format, hg19. Canonical variant input for integration tests.
- **`variants_hg38.tsv`** — the same 500 variants lifted to hg38. Used for liftover integration tests (hg38 variant list against hg19 VCFs).
- **`traits.tsv`** — 5 traits in `trait_id | trait_name | N | K | vcf_path | build` format (header row required). `ieu-a-7` (Coronary artery disease) is the binary trait and requires a `K` value; the four continuous traits leave `K` empty.
- **`vcf/`** — 5 GWAS-VCF files (hg19, CSI-indexed), one per trait.

## Flagged ambiguities

- **"REF" is overloaded**: in VCF files it means the genome-reference allele; in pleiodb it means A1 (the alphabetically-first allele). When discussing a VCF field, say "genome-ref allele"; when discussing pleiodb storage say "A1" or "A2" to avoid confusion.
- **Symmetry of p(1−p) does not imply symmetry of effect direction**: the var_y formula uses `p(1−p)` which is symmetric, but z-scores and betas are signed with respect to the A2 (ALT/effect) allele count. Saying "it doesn't matter which allele frequency you use" applies only to variance estimation — never to effect size sign.
