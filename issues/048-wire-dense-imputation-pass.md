## Parent PRD

`issues/prd-dense-vcf-imputation.md`

## What to build

Wire the dense imputation worker into the existing post-build imputation pass in `build.py`. After the V×T matrix is written, the imputation pass reads VCF paths from the database's `traits.tsv`, passes them to `impute_z_block` (which dispatches to the dense block workers), and rewrites the matrix atomically. The build phase itself is unchanged — this is purely a change to what the post-build pass hands to the workers.

See "Imputation phase — block-level design" and the overall solution description in the parent PRD.

## Acceptance criteria

- [ ] `impute_z_block` accepts `vcf_paths: list[tuple[str, str|None]]` (vcf_path, vcf_build per trait) and `vcf_threads: int`; passes them through to each block worker
- [ ] The post-build pass in `build_database` reads `vcf_path` and `vcf_build` from the just-written `traits.tsv` and supplies them to `impute_z_block`
- [ ] Blocks where a VCF read fails log a warning and fall back to sparse-only imputation for that block (no crash)
- [ ] The atomic matrix rewrite (zscore_tmp → zscore, neff_tmp → neff) is preserved unchanged
- [ ] Integration test: 5-trait build with toy LD panel and real VCF fixtures; assert imputed cell count > 0, imputed z-scores finite, observed z-scores unchanged vs reference build
- [ ] Integration test: assert imputed z-scores from dense pass are closer (higher Pearson r) to reference values than sparse imputation would produce on the same data

## Blocked by

- Blocked by `issues/044-vcf-paths-in-traits-tsv.md`
- Blocked by `issues/047-dense-imputation-worker.md`

## User stories addressed

- User story 1
- User story 2
- User story 3
- User story 5
- User story 6
- User story 9
