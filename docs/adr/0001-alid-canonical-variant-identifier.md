# ALID as the canonical variant identifier

Variants are identified by `CHROM:POS_A1_A2` (ALID) where A1 < A2 alphabetically, rather than by rsID or genome-reference-aligned `CHROM:POS:REF:ALT`. This encoding is reference-genome-free: all coordinates (chromosome, position, both alleles) are parseable from the string alone, and the canonical orientation is determined by alphabetical order of the alleles rather than by which allele matches the genome reference sequence.

## Considered options

- **rsID** — requires annotation lookup; not universally available for all variants.
- **`CHROM:POS:REF:ALT` (genome-ref aligned)** — requires a reference genome to define REF; orientation varies across GWAS studies and builds.
- **ALID (alphabetical)** — reference-free, unambiguous, computable from any two-allele record regardless of source. Chosen.

## Consequences

When ingesting a GWAS-VCF whose genome-ref REF > ALT alphabetically, the z-score is negated on read. The effect allele frequency (EAF) stored is always the frequency of A2 (the alphabetically-second allele), which is the consistent effect allele across all traits.
