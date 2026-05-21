# EAF (effect allele frequency) replaces RAF

The per-variant allele frequency stored in pleiodb is the frequency of the effect allele (A2 — the alphabetically-second allele), not the reference allele frequency (RAF) that most GWAS tools conventionally report. This follows directly from the ALID convention (ADR-0001): since A2 is always the effect allele, its frequency is the unambiguous quantity for beta/SE reconstruction.

## Consequences

The `--raf` CLI flag is removed. The variant input file (`--variants`) is a two-column TSV of `ALID | EAF`. When ingesting a VCF in non-canonical orientation (genome-ref REF > ALT alphabetically), EAF = 1 − VCF `AF`. The `raf.f16` binary file in the database is renamed to `eaf.f16`.
