## Parent PRD

`issues/prd-pleiodbr.md`

## What to build

Internal format-reading functions (not exported) that all query functions
build on:

- `read_cidx(path)` — read uint64 offset array from a `.cidx` file using
  `readBin(..., size=8, signed=FALSE)`; return as numeric (safe up to 9 PB)
- `read_chunk_raw(bin_path, cidx, chunk_id)` — seek to offset, read raw bytes
- `decompress_chunk(raw_bytes)` — zstd decompress via CRAN `zstd` package
- `decode_z(int16_vec)` — divide by 100, replace −32768 with `NA_real_`
- `decode_neff(uint16_vec)` — log2-decode back to Neff float
  (`2^(uint16 / 2048)`, NA when value is `0xFFFF`)
- `chunk_id(vi, ti, n_v_chunks)` — `ti * n_v_chunks + vi`
- `chunk_bounds(vi, ti, CV, CT, V, T)` — returns `list(v0, v1, t0, t1)`
- `get_block(db, matrix_name, v_start, v_end, t_start, t_end)` — reads all
  overlapping chunks and assembles a decoded numeric matrix

## Acceptance criteria

- [ ] `get_block(db, "zscore", 0, 512, 0, 512)` returns a 512×512 numeric
      matrix with NAs where sentinel values were
- [ ] `get_block` works when the requested range spans multiple chunks
- [ ] `decode_neff` round-trips within 0.1% of the stored value
- [ ] All functions are unexported (`.` prefix or in `R/internal.R`)

## Blocked by

- `issues/051-pleiodbr-package-scaffold.md`

## User stories addressed

None directly — internal foundation for all query functions.
