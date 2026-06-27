## Parent PRD

`issues/prd-dense-vcf-imputation.md`

## What to build

New function in `impute.py` that reads a block's `.tsv` file (already used by `build_block_index`) and returns a lookup from allele ID to LD matrix row index. Both allele orientations are included so that allele-flipped VCF variants can be matched. This lookup is used by the dense imputation worker to map z-scores returned by `read_vcf_region` to their row in the LD eigenvector matrix.

See "LD panel variant matching within block" in the parent PRD.

## Acceptance criteria

- [ ] `build_block_allele_lookup(tsv_path) -> dict[str, tuple[int, bool]]` returns `{allele_id: (ld_row_index, is_flipped)}` for every variant in the TSV
- [ ] Both orientations are present: `CHR:POS_OA_EA` → `(i, False)` and `CHR:POS_EA_OA` → `(i, True)`
- [ ] When `is_flipped=True`, callers must negate the z-score before using it as a training value
- [ ] Unit test: read a real block TSV from the LD reference panel fixture and assert correct row indices and flip flags for a handful of known variants

## Blocked by

None — can start immediately.

## User stories addressed

- User story 8
- User story 10
