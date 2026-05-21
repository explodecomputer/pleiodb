# pleiodb — Storage and Query Specification

## Overview

`pleiodb` stores GWAS summary statistics for **V = 100,000 causal variants**
across **T = 20,000 traits** in a compact, chunk-compressed binary format.
The design prioritises:

1. **Minimum on-disk footprint** — quantisation + Zstandard compression
2. **Fast random access** — pre-chunked layout, offset index for O(1) chunk seek
3. **Full reconstructibility** — beta and SE can be recovered from stored quantities
4. **Query flexibility** — variant range, trait column, p-value, arbitrary block

---

## Directory layout

```
study.pleiodb/
├── meta.json           ← dimensions, chunk sizes, encoding params, format version
├── variants.npy        ← structured array (V rows)
├── traits.npy          ← structured array (T rows)
├── zscore.bin          ← V×T  int16  compressed chunks
├── zscore.cidx         ← uint64 offset index for zscore.bin
├── neff.bin            ← V×T  uint16 compressed chunks
├── neff.cidx
├── eaf.f16             ← V    float16  effect allele freq (A2, one value per variant)
├── neff_base.f32       ← T    float32  (per-trait median Neff, fallback)
├── lambda.bin          ← T×T  float16  compressed chunks
├── lambda.cidx
└── masks/
    ├── 5e-8.coo.zst    ← COO pairs (v_idx u32, t_idx u32), Zstd-compressed
    └── 1e-5.coo.zst
```

---

## Matrices and encodings

### 1. Z-score matrix  `zscore`  — V × T, `int16`

| Property | Value |
|---|---|
| Encoding | `stored = round(z × 100)`, clamped to `[−32767, 32767]` |
| NA sentinel | `−32768` (INT16\_MIN) |
| Precision | 0.01 in z-score |
| Effective range | ±327 (covers every plausible GWAS z-score) |
| Raw size | V×T×2 = **4 GB** |
| Typical compressed | **~2.5–3.5 GB** |

**Rationale:** int16 at scale 100 gives two decimal places of z-score precision.
This translates to sub-percent error in p-values for all practically relevant
effect sizes.  int8 was rejected because it cannot faithfully represent
z-scores above ~10 with meaningful precision.

**Reconstruction:**
```
z = stored / 100       (NA where stored == −32768)
```

---

### 2. Effective sample size  `neff`  — V × T, `uint16`

| Property | Value |
|---|---|
| Encoding | `stored = round(log2(Neff) × 2048)` |
| NA sentinel | `0xFFFF` (65535) |
| Precision | ~0.034 % relative error (11 fractional bits of log2) |
| Range | 2^0 to 2^31 (~4 billion) |
| Raw size | V×T×2 = **4 GB** |
| Typical compressed | **~1.5–2.5 GB** (Neff is correlated within traits) |

**Rationale:** logarithmic uint16 encoding gives better relative precision than
float16 across the full range of sample sizes (hundreds to millions).  The
0.034 % relative error in Neff propagates to ~0.017 % error in SE, which is
negligible for any downstream application.

**Reconstruction:**
```
Neff = 2^(stored / 2048)     (NA where stored == 0xFFFF)
```

---

### 3. Effect allele frequency  `eaf`  — V × 1, `float16`

Frequency of A2 (the effect / alphabetically-second allele).  Stored as a
flat binary file of `V` float16 values (200 KB raw).  float16 gives ~0.1 %
relative precision for frequencies in (0.001, 0.999), sufficient for SE
reconstruction.

---

### 4. Sample overlap / lambda matrix  `lambda`  — T × T, `float16`

Estimated as pairwise Pearson correlations of z-scores at null variants
(|z| < 3.0 in both traits).  Stored symmetric (full T×T) in the same
chunked format as the main matrices.

| Property | Value |
|---|---|
| Shape | T×T = 400 M cells |
| dtype | float16 |
| Raw size | **800 MB** |
| Typical compressed | **~400–600 MB** |

float16 is sufficient for correlation values (range −1 to 1, precision ~0.001).

---

### 5. Significance masks  `masks/`  — sparse, per threshold

Each mask is a Zstd-compressed binary file containing a flat array of
`(v_idx: uint32, t_idx: uint32)` pairs sorted lexicographically by
`(v_idx, t_idx)`.

| Threshold | Expected hits | Size |
|---|---|---|
| 5×10⁻⁸ (GWS) | ~0.1–1 M | **0.6–6 MB** compressed |
| 1×10⁻⁵ | ~5–50 M | **30–300 MB** compressed |

Masks allow O(hits) p-value queries without scanning the full matrix.
When a mask is not present for a requested threshold, the query falls back
to a full scan.

---

## Chunk format

Every chunked matrix (zscore, neff, lambda) uses an identical two-file layout:

```
{name}.bin   — concatenated Zstandard-compressed chunk blobs
{name}.cidx  — uint64 array, length = n_chunks + 1
               cidx[i]   = start byte of chunk i in .bin
               cidx[i+1] = end byte of chunk i in .bin
```

**Chunk ordering:** row-major over (v\_chunk, t\_chunk).
```
chunk_id = vi * n_t_chunks + ti
```

**Default chunk dimensions:** 512 × 512 (configurable at build time).

**Reading chunk (vi, ti):**
1. `start = cidx[vi * n_t_chunks + ti]`
2. `end   = cidx[vi * n_t_chunks + ti + 1]`
3. `seek(start); read(end − start)` → Zstd decompress → reshape to int16/uint16

Chunks at the boundary of the matrix are smaller (clipped to actual dimensions).

---

## Compression

**Algorithm:** Zstandard (zstd), level 3.
Chosen over LZ4 for better ratio with similar decompression speed.
Blosc was considered but rejected to avoid the HDF5/Blosc dependency chain.

**Expected compression ratios:**

| Matrix | Ratio | Notes |
|---|---|---|
| zscore int16 | 1.2–1.5× | Near-random; shuffle helps marginally |
| neff uint16 | 2–4× | Many identical values within a trait |
| lambda float16 | 1.5–2× | Smooth correlation structure |
| masks COO | 3–10× | Very sparse at GWS threshold |

---

## Storage budget (V=100k, T=20k)

| Component | Raw | Compressed (est.) |
|---|---|---|
| Z-score | 4.0 GB | 2.5–3.0 GB |
| Neff | 4.0 GB | 1.5–2.0 GB |
| RAF | 0.0002 GB | negligible |
| Lambda | 0.8 GB | 0.4–0.6 GB |
| Sig masks (2 thresh) | 0.05–0.6 GB | 0.03–0.3 GB |
| Metadata / indices | <10 MB | — |
| **Total** | **~9 GB** | **~5–6 GB** |

---

## Input file formats

### `--variants`  (required)

Tab-separated, no header.  Lines beginning with `#` are ignored.

| Column | Type | Description |
|--------|------|-------------|
| 1 `alid` | string | **ALID** — `CHROM:POS_A1_A2` where A1 ≤ A2 alphabetically (e.g. `10:101558746_G_T`). Encodes all coordinates; must be unique. |
| 2 `eaf` | float | Effect allele frequency for A2 (0 < EAF < 1). Optional: omit or leave blank to store NaN (beta/SE reconstruction will not be possible). |

Variants with non-canonical allele order (A1 > A2) are automatically
normalised on load: the ALID is rewritten and EAF is flipped to `1 − EAF`.

VCF records are matched by **position** (CHROM:POS), then confirmed by
allele pair.  When the VCF alleles are non-canonical (genome-ref REF > ALT
alphabetically), the z-score is negated so it reflects the effect of A2.

**Example:**
```
10:101558746_G_T	0.56079
10:1206798_C_T	0.071479
19:44908822_C_T	0.12
```

---

### `--traits`  (required)

Tab-separated, no header.  Lines beginning with `#` are ignored.

| Column | Type | Description |
|--------|------|-------------|
| 1 `trait_id` | string | Unique identifier for the trait (e.g. `ieu-b-2`, `ukb-b-1234`). Stored in the database and used in all query output. |
| 2 `trait_name` | string | Optional human-readable label (e.g. `Body mass index`). Stored in the database; empty string if absent. |
| 3 `vcf_path` | string | Absolute or relative path to a GWAS-VCF file (bgzipped `.vcf.gz` with CSI or TBI index recommended). |
| 4 `build` | string | Optional genome build of the VCF (`hg19`, `hg38`, `GRCh37`, `GRCh38`). Used with `--variants-build` to trigger automatic liftover when builds differ. |

**Example:**
```
ieu-b-2	Body mass index	/data/gwas/body_mass_index.gwas.vcf.gz	hg19
ukb-b-1234	LDL cholesterol	/data/gwas/ukbb_ldl.gwas.vcf.gz	hg38
finn-r-T2D		/data/gwas/finngen_T2D.vcf.gz
```

---

### GWAS-VCF files  (one per trait)

pleiodb reads the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification)
(Lyon et al. 2021).  The following FORMAT fields are consumed:

| Field | Required | Description |
|-------|----------|-------------|
| `ES`  | Yes* | Effect size (beta). Used to compute z = ES / SE. |
| `SE`  | Yes* | Standard error. |
| `EZ`  | No | Pre-computed z-score. Used instead of ES/SE when present. |
| `SS`  | No | Per-variant sample size. Stored as Neff; falls back to `NS` INFO field. |

\* Either `EZ` or both `ES` + `SE` must be present.

Variants in the VCF are matched to the variant list by the VCF `ID` column
(default) or by `CHROM:POS:REF:ALT` key (pass `--id-col CHRPOSREFALT` in
future versions).  Variants in the VCF that do not appear in the variant list
are silently skipped.  Variants in the list that are absent from a VCF receive
`NaN` for that trait.

---

## Beta / SE reconstruction

Given z-score, Neff, and EAF for variant v in trait t:

```
SE   = 1 / sqrt(Neff × 2 × p × (1 − p))
beta = z × SE
```

where `p = EAF` (effect allele frequency of A2).  This formula is symmetric
in p and 1−p, so EAF and RAF give identical SE; only the sign of beta differs
depending on which allele is the effect allele — in pleiodb it is always A2.

Reconstruction error budget:
- z quantisation: ±0.005 → ±0.5 % error in beta at z=1
- Neff quantisation: ±0.034 % → negligible SE error
- EAF (float16): ±0.1 % → negligible SE error

---

## GWAS-VCF ingestion

Supports the GWAS-VCF specification (Lyon et al. 2021, *Nature Genetics*).

Required FORMAT fields:
- `ES` — effect size (beta)
- `SE` — standard error

Optional FORMAT fields (used when present):
- `EZ` — pre-computed z-score (preferred over ES/SE ratio)
- `SS` — per-variant sample size

**Variant matching:** by position (CHROM:POS), then allele pair (REF, ALT).
The VCF `ID` field (rsID) is not used.  When the VCF has genome-ref REF > ALT
alphabetically (non-canonical orientation), the z-score is negated so that the
stored value always reflects the effect of A2 (the alphabetically-second allele).

**bcftools pre-filter:** when the VCF has a CSI or TBI index, `bcftools view -R`
is used to restrict the file to queried positions before parsing.  Falls back
to a full cyvcf2 scan if bcftools is unavailable or the index is absent.

---

## Query patterns and complexity

| Query type | Chunks read | Notes |
|---|---|---|
| Single variant, all traits | `n_t_chunks = ceil(T/512)` | 40 chunks |
| Single trait, all variants | `n_v_chunks = ceil(V/512)` | 196 chunks |
| Genomic region (k variants) | `ceil(k/512) × n_t_chunks` | sorted by offset |
| p-value filter | O(mask_hits) | uses pre-built COO mask |
| V × T block | `ceil(|V|/512) × ceil(|T|/512)` | arbitrary sub-matrix |

Chunks within a query are sorted by their byte offset before reading to
maximise sequential I/O throughput on spinning disk and SSD alike.

---

## Build process

```
pleiodb build  OUTPUT_DIR \
  --variants variants.tsv \   # ALID  EAF  (CHROM:POS_A1_A2, alphabetical alleles)
  --traits   traits.tsv   \   # trait_id  trait_name  vcf_path  [build]
  --workers  16            \
  --chunk-v  512  --chunk-t 512
```

**Algorithm:**
1. Load variant list → integer index (variant → row).
2. Iterate traits in batches of `chunk_t` (one "column-slab").
3. Within each batch, read VCF files in parallel using a thread pool.
4. Quantise the `V × batch_T` z-score and Neff sub-matrices.
5. Write row-chunks (size `chunk_v × batch_T`) to the .bin/.cidx files.
6. After all traits: write RAF, compute and write significance masks.
7. Optionally compute lambda matrix with `pleiodb lambda`.

Peak memory per process: `V × batch_T × 8 bytes ≈ 400 MB` at default settings.

---

## Lambda computation

```
pleiodb lambda  STUDY.pleiodb  --null-thresh 3.0  --workers 16
```

For each pair of traits (t1, t2), Pearson correlation is computed over the
set of variants where both |z\_t1| < 3.0 and |z\_t2| < 3.0.  A minimum of
30 variants is required; otherwise the entry is set to NaN.

The lambda matrix is written as a separate chunked float16 matrix.

---

## CLI reference

```
pleiodb build   OUTPUT  --variants TSV  --traits TSV  [options]
pleiodb lambda  DB      [--null-thresh FLOAT]  [--workers N]
pleiodb query   DB      [query options]  [--format tsv|json]  [-o FILE]
pleiodb info    DB
```

### Query options

| Flag | Description |
|---|---|
| `--variant ID` | All T z-scores for one variant |
| `--trait ID` | All V z-scores for one trait |
| `--region chr:start-end` | Variants in region × all traits |
| `--variants-file FILE` | Subset of variants |
| `--traits-file FILE` | Subset of traits |
| `--pval FLOAT` | Filter by p-value (uses mask if available) |
| `--beta-se` | Add reconstructed beta and SE columns |
| `--format tsv\|json` | Output format (default: tsv) |

---

## Python API

```python
import pleiodb

db = pleiodb.GWASDatabase.open("study.pleiodb")

# All traits for one variant
z = db.zscore_variant(v_idx=0)                         # float32 (T,)

# All variants for one trait
z = db.zscore_trait(t_idx=5)                           # float32 (V,)

# Genomic region
v_idx, z_mat = db.zscore_region("chr1", 1_000_000, 2_000_000)  # (k,), (k,T)

# Arbitrary block
z = db.zscore_block(v_indices=[0,1,2], t_indices=[10,11])       # (3, 2)

# Reconstruct beta and SE
beta, se = db.beta_se_variant(v_idx=0)                          # float32 (T,), (T,)
beta, se = db.beta_se_block([0,1], [10,11])                     # (2,2), (2,2)

# Significant associations
v_idx, t_idx, z_vals = db.query_significant(pval=5e-8)
```

---

## Design decisions and trade-offs

### Why not HDF5?

HDF5 was considered and rejected for the following reasons:
- Requires `libhdf5` system library, which creates deployment friction.
- SWMR (single-writer multiple-reader) mode is complex to configure correctly.
- The HDF5 chunking model is equivalent to what this format implements directly,
  without the additional abstraction layer.

The custom format achieves the same chunking + compression performance with
~200 lines of code and a single Python dependency (`zstandard`).

### Why not Zarr?

Zarr (v2) creates one file per chunk, leading to ~16,000 files for V×T at
default chunk sizes — problematic on HPC filesystems with inode quotas.
Zarr v3 supports "sharding" but is still experimental.

### Why int16 and not int8 for z-scores?

int8 with a clip at ±10 gives 0.08 precision in z-score.  This translates to
~4 % error in p-value near the genome-wide significance threshold (z ≈ 5.45).
For applications like Mendelian randomisation, colocalization, and PheWAS,
this error is unacceptable.  int16 at scale 100 gives 0.01 precision with
only 2× the storage cost.

### Why uint16 log2-encoding for Neff?

Linear uint16 would limit Neff to 65,535, excluding large biobank studies.
float16 was considered but has non-uniform precision (coarser at large N).
Log2 encoding with 11 fractional bits gives 0.034 % relative error uniformly
across the range 1–2×10⁹, covers all plausible sample sizes, and fits in uint16.

### Why store the full T×T lambda matrix (not upper triangle)?

Upper-triangle storage saves 50 % space (~400 MB) but requires branch logic
on every access and complicates the chunk-based random-access pattern.
At ~400–600 MB compressed, the full symmetric storage is acceptable.

### Chunk size choice (512 × 512)

- Single-variant query: 40 chunks (~8 MB compressed data) — fast
- Single-trait query: 196 chunks (~40 MB compressed data) — acceptable
- Chunk file size: ~200 KB after compression — large enough for sequential I/O
  efficiency, small enough for granular random access

A 256×1000 chunk would favour trait queries; 1000×256 would favour variant
queries.  512×512 is a symmetric compromise.  The chunk size is configurable
at build time, so users with skewed query patterns can optimise accordingly.

### Significance masks

Pre-built COO masks enable O(hits) p-value queries without scanning 2 billion
cells.  Two thresholds (5×10⁻⁸ and 10⁻⁵) cover genome-wide significance and
a permissive suggestive threshold.  Additional thresholds can be added without
rebuilding the main matrices.
