## Parent PRD

`issues/prd-dense-vcf-imputation.md`

## What to build

Modify `_impute_block_process` in `impute.py` to accept dense z-scores from the full LD panel as the elastic net training set, while keeping the missing pleiodb positions as the only prediction targets. The worker also gains the ability to read VCF regions in parallel using an inner ThreadPoolExecutor: given a list of (trait_index, vcf_path) pairs for traits with missing data in this block, it reads all of their VCF regions concurrently (threads, since bcftools is subprocess-based and not GIL-limited), assembles the dense z-score matrix, then fits and predicts as before.

See "Imputation phase — block-level design" and "Parallelism — nested two-level" in the parent PRD.

## Acceptance criteria

- [ ] Worker accepts `vcf_paths: list[tuple[int, str]]` (trait_index, vcf_path) and `vcf_threads: int` in addition to existing args
- [ ] VCF region reads for missing traits are dispatched via `ThreadPoolExecutor(max_workers=vcf_threads)` within the worker process
- [ ] Only traits with at least one NaN in `z_pleiodb[:, j]` for this block issue a VCF read; fully-observed traits are skipped
- [ ] Dense z-scores (from VCF, covering all LD panel variants in the block) are used as the elastic net training set; observed pleiodb z-scores are NOT used for training
- [ ] Predictions are written only to missing pleiodb positions (no overwrite of observed data)
- [ ] When `vcf_paths` is empty or None, behaviour is identical to the current sparse worker (backward-compatible)
- [ ] Unit test: synthetic block where pleiodb covers 5/100 LD panel variants; dense training (100 variants) produces Pearson r > 0.8 vs held-out truth; sparse training (5 variants) produces r < 0.8
- [ ] Unit test: assert observed pleiodb z-scores are identical before and after the dense worker runs

## Blocked by

- Blocked by `issues/045-read-vcf-region-function.md`
- Blocked by `issues/046-block-tsv-allele-lookup.md`

## User stories addressed

- User story 1
- User story 3
- User story 5
- User story 7
- User story 8
- User story 10
