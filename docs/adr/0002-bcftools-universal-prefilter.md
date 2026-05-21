# bcftools used as the universal VCF pre-filter

Every GWAS-VCF read uses `bcftools view -R <positions>` to restrict the file to the queried positions before cyvcf2 parses individual records, not only in the liftover path. bcftools is therefore a required runtime dependency (installed via conda/bioconda).

## Considered options

- **cyvcf2 full scan** — simple, no extra dependency. Unacceptably slow for production VCFs with millions of variants; fine only for the small test files.
- **tabix/cyvcf2 region queries** — fast but requires one query call per variant; overhead dominates for 100k variants.
- **bcftools pre-filter (chosen)** — one subprocess call produces a position-restricted VCF; cyvcf2 then does a fast linear scan of the small result. Symmetric behaviour for liftover and non-liftover paths. Falls back to full cyvcf2 scan if no index is present or bcftools is unavailable.
