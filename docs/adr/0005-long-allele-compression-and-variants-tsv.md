# Compress long alleles in ALIDs; store variants as TSV not binary

## Variant storage

`variants.npy` (structured array with U512 allele fields) and `eaf.f16` are replaced by a human-readable `variants.tsv` inside the `.pleiodb` directory with columns `ALID` and `EAF`. At 1M variants with typical SNP alleles (~1 char each), the numpy U512 dtype allocates 2 048 bytes per allele — ~4.3 GB total — while the TSV is ~30 MB. `chrom`, `pos`, `a1`, `a2` are derived by parsing the ALID on load and are not stored separately.

## Long-allele compression

Alleles longer than 20 characters are compressed to `{allele[:8]}~{sha256(allele)[:4]}` (13 chars) in both the ALID and the `variants.tsv`. The tilde (`~`) is the truncation marker — not a valid nucleotide character. Four hex chars of SHA-256 (16 bits) are sufficient for collision resistance because uniqueness is only required among alleles at the same genomic position, where at most a handful of distinct long alleles ever appear.

## Consequences

- ALIDs for long-allele variants are not fully reversible to the original sequence without the source VCF.
- Any tool that generates ALIDs must apply the same compression rule for alleles > 20 chars.
