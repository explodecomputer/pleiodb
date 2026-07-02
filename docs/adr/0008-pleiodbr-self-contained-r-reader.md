# ADR 0008 — pleiodbr: self-contained R reader for .pleiodb datasets

**Status**: accepted

## Context

R users need to query `.pleiodb` datasets without installing Python or the
`pleiodb` Python package. Two approaches were considered:

1. **reticulate wrapper** — thin R shim that delegates to the Python `GWASDatabase`
   class via `reticulate`.
2. **Native R reader** — re-implement the format parsing directly in R.

## Decision

Build `pleiodbr` as a **self-contained native R package** that reads `.pleiodb`
files directly, with no Python dependency.

The binary format is simple enough to implement in R:

- `meta.json` — standard JSON
- `variants.tsv` / `traits.tsv` — standard TSV
- `zscore.bin` / `zscore.cidx` — zstd-compressed int16 chunks with uint64 offset
  index; decompressed via the CRAN `zstd` package
- `imputed.coo.zst` — zstd-compressed uint32 pair array

## API

All query functions take an S3 connection object as their first argument
(functional style, compatible with `|>`).

```r
db <- open_pleiodb("/path/to/main.pleiodb")

phewas(db, "1:103574777_C_G")      # single ALID
phewas(db, "1:103e6-104e6")        # region string chrom:start-end (bp)
gwas(db, "ukb-b-19953")
tophits(db, traits, pval = 5e-8)   # traits required
associations(db, variants, traits)
rho(db, traits_1, traits_2)
```

All functions return a **tibble** with a consistent column schema:

| Column | Type | Notes |
|--------|------|-------|
| `variant_id` | character | ALID (`chrom:pos_A1_A2`) |
| `trait_id` | character | OpenGWAS ID |
| `z` | double | |
| `beta` | double | reconstructed from z + neff + eaf |
| `se` | double | |
| `pval` | double | |
| `eaf` | double | effect allele frequency |
| `n` | double | effective N |
| `imputed` | logical | from `imputed.coo.zst` |

`rho()` returns a tibble with columns `trait_id_1`, `trait_id_2`, `rho` (long
format; use `tidyr::pivot_wider()` for matrix form).

## Consequences

- No Python dependency; installs cleanly from CRAN or GitHub.
- Format parsing must be kept in sync with the Python implementation if the
  `.pleiodb` format changes. The format version in `meta.json` should be checked
  on open and an informative error raised if it is newer than the package
  supports.
- `reticulate` wrapper is not built; R users who want Python-only features (build,
  imputation, rho computation) must use the Python package directly.
