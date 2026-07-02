# PRD — pleiodbr: self-contained R reader for .pleiodb datasets

## Background

`pleiodb` is a chunk-compressed binary store for GWAS z-scores (V×T matrix,
95,378 variants × 4,159 traits). The format is defined by the Python package.
R users need to query `.pleiodb` datasets without installing Python.

See also: `docs/adr/0008-pleiodbr-self-contained-r-reader.md`

## Goal

Ship `pleiodbr`, a standalone R package (own GitHub repo,
`github.com/explodecomputer/pleiodbr`) that reads `.pleiodb` directories
natively. No Python or `reticulate` dependency.

## User stories

1. As an R user, I can open a `.pleiodb` dataset and inspect its dimensions and
   metadata without writing any boilerplate.
2. As an R user, I can run a PheWAS by ALID string and get a tidy tibble of
   z-scores, betas, SEs, p-values, EAFs, Neff, and imputation status across
   all traits.
3. As an R user, I can run a PheWAS for all variants in a genomic region by
   passing a `"chrom:start-end"` string.
4. As an R user, I can pull a full GWAS (all variants for one trait) as a
   tidy tibble.
5. As an R user, I can retrieve the top hits (p < threshold) for one or more
   specified traits as a tidy tibble.
6. As an R user, I can query an arbitrary block of variants × traits as a tidy
   tibble.
7. As an R user, I can retrieve pairwise rho (phenotypic correlation) values
   for sets of traits as a tidy tibble.
8. As an R user, I get an informative error if the `.pleiodb` format version is
   newer than the package supports.

## API contract

### Connection

```r
db <- open_pleiodb("/path/to/main.pleiodb")
# S3 object of class "pleiodb"
# print(db) shows V, T, format version, path
```

### Query functions

All return a **tibble** unless noted.

```r
# User stories 2 & 3
phewas(db, variant)          # variant = ALID string OR "chrom:start-end"

# User story 4
gwas(db, trait_id)

# User story 5
tophits(db, traits, pval = 5e-8)   # traits required

# User story 6
associations(db, variants, traits)  # vectors of ALIDs and trait IDs

# User story 7
rho(db, traits_1, traits_2)        # long format: trait_id_1, trait_id_2, rho
```

### Tibble schema (all functions except rho)

| Column | Type | Notes |
|---|---|---|
| `variant_id` | character | ALID (`chrom:pos_A1_A2`) |
| `trait_id` | character | OpenGWAS ID |
| `z` | double | |
| `beta` | double | reconstructed: `z / sqrt(neff * 2*eaf*(1-eaf))` |
| `se` | double | `1 / sqrt(neff * 2*eaf*(1-eaf))` |
| `pval` | double | two-sided normal |
| `eaf` | double | effect allele frequency |
| `n` | double | effective N |
| `imputed` | logical | TRUE if from elastic-net imputation pass |

### rho tibble schema

| Column | Type |
|---|---|
| `trait_id_1` | character |
| `trait_id_2` | character |
| `rho` | double |

## Format details (for implementation)

| File | Content |
|---|---|
| `meta.json` | V, T, chunk shape, format_version, pval_thresholds |
| `variants.tsv` | alid, eaf (+ optional alid_hg38) |
| `traits.tsv` | trait_id, trait_name, N, K, neff_study, … |
| `zscore.bin/.cidx` | V×T int16, z×100, NA=−32768, ti-major chunks |
| `neff.bin/.cidx` | V×T uint16, log2-encoded, NA=0xFFFF |
| `rho.bin/.cidx` | T×T float16, same chunk layout |
| `imputed.coo.zst` | zstd-compressed uint32 pairs (v_idx, t_idx) |
| `masks/{thr}.coo.zst` | significance hit pairs for `tophits()` |

Chunk layout: chunk_id = `ti * n_v_chunks + vi`. cidx is a uint64 array of
length `n_chunks + 1`; `cidx[i]` and `cidx[i+1]` give the byte range in `.bin`.

## Dependencies

- `zstd` (CRAN) — chunk decompression
- `tibble` — return type
- `dplyr` — internal data manipulation

## Out of scope (v1)

- Build / imputation / rho computation (Python-only)
- rsID lookup
- CRAN submission (GitHub-only for now)
