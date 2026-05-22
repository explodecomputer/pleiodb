# Derive Neff from SE and traits.tsv; ignore VCF SS field

Per-variant Neff is computed at ingest as `var_y[t] / (SE[v,t]² × 2·EAF[v]·(1−EAF[v]))` rather than being read from the VCF `SS` field. `var_y[t]` is itself estimated as the median of `SE[v,t]² × 2·EAF[v]·(1−EAF[v]) × Neff_study[t]` across common variants (EAF ∈ (0.01, 0.99)), where `Neff_study[t]` = `N[t]` for continuous traits and `4·N[t]·K[t]·(1−K[t])` for binary traits, with `N` and `K` (case fraction) supplied in the traits input file.

## Why not use VCF SS

The `SS` FORMAT field is not consistently defined across GWAS software: some tools write total sample size (`Ncase + Ncontrol`), some write binary effective N (`4/(1/Ncase + 1/Ncontrol)`), and some omit it entirely. Using raw SS as Neff in the normalised-beta reconstruction formula (`SE_norm = 1/sqrt(Neff × 2p(1−p))`) requires that SS equals `N / var_y` — an assumption that holds only incidentally and cannot be verified from the VCF alone.

## Considered alternatives

- **Read SS from VCF directly** — rejected because SS semantics are inconsistent across tools (see above).
- **Require users to supply Neff explicitly** — rejected because Neff for binary traits depends on case fraction in a non-obvious way; N + K is a more natural and auditable input.

## Consequences

- `ES` and `SE` FORMAT fields are required in every ingested VCF (used to derive both Z and Neff).
- `N` is a required column in the traits input file; `K` is required for binary traits.
- The stored Neff is on the normalised (var(y)=1) scale: it equals the actual study N only when the GWAS formula holds exactly and EAF_user ≈ EAF_study.
- The single shared EAF per variant (from the user variants file) introduces a small approximation in Neff and var_y when study-specific EAFs differ from it.
