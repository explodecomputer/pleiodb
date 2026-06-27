## Parent PRD

`issues/prd-dense-vcf-imputation.md`

## What to build

New function `read_vcf_region(vcf_path, chrom, start, end)` in `vcf.py` that extracts z-scores for **all** variants in a genomic region from a GWAS-VCF file. Returns a dict `{allele_id: z_score}` where allele_id matches the format used in the LD panel TSV (`CHR:POS_OA_EA`). Unlike the existing `read_vcf`, this function is not filtered to the 95k pleiodb variant list — it returns everything bcftools finds in the region.

See "VCF region reading — new function" in the parent PRD.

## Acceptance criteria

- [ ] `read_vcf_region(vcf_path, chrom, start, end)` returns a `dict[str, float]` keyed by allele ID for all variants found in the region
- [ ] Uses `bcftools view -r {chrom}:{start}-{end}` when the VCF has a CSI or TBI index; falls back to full cyvcf2 scan otherwise
- [ ] Z-score extracted from FORMAT/EZ if present, else FORMAT/ES / FORMAT/SE (same logic as existing `_extract_z`)
- [ ] Allele flip is handled: if VCF REF/ALT are reversed relative to the allele ID, z-score sign is flipped
- [ ] Returns an empty dict (not an error) when the region contains no matching variants
- [ ] Unit test with an existing fixture VCF: query a known region and assert the returned dict contains the expected variants and z-scores
- [ ] Unit test: querying a region with no variants returns `{}`

## Blocked by

None — can start immediately.

## User stories addressed

- User story 7
- User story 8
- User story 9
