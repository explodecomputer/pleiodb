## Parent PRD

`issues/prd.md`

## What to build

Implement `build_block_index()` in `src/pleiodb/impute.py`. This function scans the LD
reference panel directory tree once at the start of a build, matches pleiodb variant ALIDs to
the SNP IDs in each block's `.tsv` metadata file, and returns a compact mapping that the
imputation loop can reuse across all trait batches without re-scanning the filesystem.

**Directory layout expected:**
```
ld_dir/
  {chr}/
    {start}-{end}/
      {start}-{end}.tsv          (columns: CHR SNP OA EA EAF BP)
      {start}-{end}.unphased.vcor1.gz   (LD matrix, no header, tab-separated)
```

**Matching logic:** LD-panel SNP IDs are already in canonical alphabetical-allele ALID format
(`chr:bp_A1_A2`), so they can be compared directly against pleiodb `variants['id']` entries
after stripping the `chr` prefix from the directory chromosome name if needed. Blocks with
fewer than 2 matched variants are excluded from the index.

**Return value:** a dict mapping each `block_dir: Path` to:
```python
{
    'variant_indices': list[int],   # indices into the pleiodb variants array
    'ld_row_indices':  list[int],   # corresponding row/col positions in the LD matrix
    'n_ld_snps':       int,         # total SNPs in the LD panel block (for subsetting)
}
```

This slice also adds the private helper `_load_ld_submatrix(block_dir, ld_row_indices)` which
reads the `.unphased.vcor1.gz` file and returns the symmetric submatrix restricted to the
given row/column indices. This helper is placed here (not in issue 033) because it performs
filesystem I/O.

## Acceptance criteria

- [ ] `build_block_index(variants, ld_dir, ancestry='EUR')` is exported from `impute.py`.
- [ ] Given a synthetic LD panel directory with two chromosome subdirectories and two blocks
      each (created by the test fixture), the function returns an index containing only the
      blocks that have ≥ 2 pleiodb variants matched to LD-panel SNP IDs.
- [ ] Blocks where zero or one pleiodb variants match the LD panel are silently excluded
      (no error or warning).
- [ ] A block directory that is missing its `.tsv` file is skipped without raising an
      exception (logged at WARNING level).
- [ ] `_load_ld_submatrix(block_dir, ld_row_indices)` returns a numpy array of shape
      `(len(ld_row_indices), len(ld_row_indices))` that equals the expected rows/cols subset
      of the full LD matrix in the fixture.
- [ ] The function completes in under 5 seconds when pointed at the real LD panel
      (`/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38/EUR`)
      with the production variant list (manual/smoke test only; not a CI assertion).
- [ ] Tests in `tests/test_impute.py` cover `build_block_index` using a small temporary
      directory fixture (created with `tmp_path`); no dependency on the real LD panel.

## Blocked by

None — can start immediately (parallel with `issues/033-imputation-kernel.md`).

## User stories addressed

- User story 2 (uses the same LD reference data as the gpmap pipeline)
- User story 6 (partial LD panel: missing blocks are skipped, not fatal)
- User story 14 (block index built once, reused across trait batches)
