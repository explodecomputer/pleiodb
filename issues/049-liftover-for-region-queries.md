## Parent PRD

`issues/prd-dense-vcf-imputation.md`

## What to build

Extend `read_vcf_region` to accept an optional `vcf_build` parameter. When `vcf_build` differs from hg38 (the LD panel coordinate system), lift the block's start/end coordinates from hg38 to the VCF's build before issuing the `bcftools view -r` call. A single interval liftover is sufficient — no per-variant lifting required.

See "Liftover for region queries" in the parent PRD.

## Acceptance criteria

- [ ] `read_vcf_region(vcf_path, chrom, start, end, vcf_build=None)` lifts `(chrom, start, end)` from hg38 to `vcf_build` when `vcf_build` is not None and not hg38
- [ ] Uses the existing liftover machinery (same chain files and `lift_variants` plumbing already in `liftover.py`)
- [ ] When liftover of the region fails (e.g. the interval spans a gap), logs a warning and returns `{}` for that block/trait combination
- [ ] When `vcf_build` is None or hg38, behaviour is identical to before (no liftover performed)
- [ ] Unit test: mock `read_vcf_region` with `vcf_build="hg19"`; assert the coordinates passed to bcftools are in hg19, not hg38

## Blocked by

- Blocked by `issues/048-wire-dense-imputation-pass.md`

## User stories addressed

- User story 6
