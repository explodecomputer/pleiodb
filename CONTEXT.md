# pleiodb

A compact binary store and query engine for GWAS summary statistics across many variants and traits simultaneously.

## Language

### Variants and alleles

**Variant**:
A biallelic locus identified by `CHROM:POS_REF_ALT` where REF and ALT are in canonical allele orientation. The identifier encodes all coordinates; no external reference genome is required.
_Avoid_: SNP (too narrow — indels are also stored), site

**Canonical allele orientation**:
The convention where REF = `min(ref, alt)` alphabetically and ALT = `max(ref, alt)`. When a source VCF has the alleles in non-canonical orientation (genome-reference REF > ALT alphabetically), the z-score is negated on ingestion to compensate.
_Avoid_: "reference-aligned", "forward strand"

**ALID (Alphabetic Variant ID)**:
The string `CHROM:POS_A1_A2` where A1 < A2 alphabetically (e.g. `10:101558746_G_T`). The canonical primary key for a variant in pleiodb. Encodes chromosome, position, and both alleles; all coordinates are parseable from the ALID alone.
_Avoid_: rsID (rsIDs are an alternate lookup key, not the ALID), "variant ID" (use ALID to be precise)

**Positional matching**:
The strategy for aligning variant IDs to VCF records: match by CHROM:POS first (fast index lookup), then confirm by REF/ALT alleles. rsIDs in the VCF `ID` field are not used for primary matching.
_Avoid_: ID-based matching

### Traits and data

**Trait**:
A GWAS phenotype represented as a column in the V×T matrix. Identified by a `trait_id` string (e.g. `ieu-a-7`).

**GWAS-VCF**:
A VCF file following the MRC IEU GWAS-VCF specification. Source of z-scores and sample sizes per variant per trait.

**Z-score**:
The primary stored statistic: effect size (ES) divided by standard error (SE), or the pre-computed `EZ` field when present.

**Neff (effective sample size)**:
Per-variant sample size stored as a log2-encoded uint16. Used for beta/SE reconstruction.

**EAF (effect allele frequency)**:
Allele frequency for the effect allele (A2 — the alphabetically second allele). The primary per-variant frequency input; stored as float16. Used for beta/SE reconstruction. When ingesting a VCF: if the VCF alleles are already in canonical orientation (A1=REF < ALT=A2), EAF = VCF `AF`; if they are flipped (REF > ALT alphabetically), EAF = 1 − VCF `AF`.
_Avoid_: RAF (retired — "reference" is ambiguous between genome-ref and canonical-ref)

### Storage

**Database** (`*.pleiodb`):
The output directory containing all binary matrices, metadata, and index files for one collection of variants × traits.

**Chunk**:
The unit of compressed storage: a rectangular submatrix of fixed size (default 512×512 cells) stored as a single Zstandard-compressed blob.

**Significance mask**:
A Zstandard-compressed COO-format sparse array of (v_idx, t_idx) pairs for all associations passing a p-value threshold. Enables O(hits) p-value queries.

## Input file formats

**Variant input file**: Two columns, tab-separated, no header: `ALID | EAF`. Provides both the variant list and per-variant effect allele frequencies in a single file.

**Trait input file**: Tab-separated, no header. Columns: `trait_id | trait_name | vcf_path | build(optional)`. `trait_name` precedes `vcf_path` for readability.

## Test data (`tests/test_data/`)

- **`variants_hg19.tsv`** — 500 variants in ALID+EAF format, hg19. Canonical variant input for integration tests.
- **`variants_hg38.tsv`** — the same 500 variants lifted to hg38. Used for liftover integration tests (hg38 variant list against hg19 VCFs).
- **`traits.tsv`** — 5 traits in `trait_id | trait_name | vcf_path` format.
- **`vcf/`** — 5 GWAS-VCF files (hg19, CSI-indexed), one per trait.

## Flagged ambiguities

- **"REF" is overloaded**: in VCF files it means the genome-reference allele; in pleiodb it means A1 (the alphabetically-first allele). When discussing a VCF field, say "genome-ref allele"; when discussing pleiodb storage say "A1" or "A2" to avoid confusion.
