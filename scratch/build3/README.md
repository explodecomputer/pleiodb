# Build 3: Full 4,159-Trait Build with Dense LD Imputation

## Goal

Full production build of the pleiodb database:
- **V = 95,378 variants** (hg19, with EAF)
- **T = 4,159 OpenGWAS phenotype traits**
- **Dense VCF imputation**: for each LD block the elastic net trains on all VCF
  variants in the block region, not just the 95k stored ones

## Inputs

| Input | Path |
|-------|------|
| Variants (hg19, with EAF) | `scratch/build1/variants_hg19_eaf.tsv` |
| Traits (all 4,159) | `scratch/build1/traits.tsv` |
| LD reference panel | `/local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38` |
| LD ancestry | `EUR` |

The `.ldeig.npz` eigenvector caches were fully populated by Build 2 so no cold
Rscript calls are needed during the imputation pass.

## Output

`/local-scratch/data/pleiodb/main.pleiodb`

## Timing Projection

Based on Build 2 (100 traits, 369 min imputation):

| Phase | Build 2 (100 traits) | Build 3 (4,159 traits) |
|-------|----------------------|------------------------|
| VCF ingestion | 2.5 min | ~115 min (9 batches × ~13 min) |
| Imputation | 369 min | ~3,200 min (35× more VCF reads/block, 4× faster with more workers) |
| **Total** | **371 min** | **~55 hours (2–3 days)** |

The dominant cost is inner VCF region reads: ~27% of 4,159 ≈ 1,123 traits need
a bcftools call per block. With 16 vcf-threads that is 70 batches vs 4 batches
for trial100 (17× more). Doubling outer workers to 32 gives ~2× speedup.

## Build Command

```bash
cd /home/gh13047/repo/pleiodb

pleiodb build \
  /local-scratch/data/pleiodb/main.pleiodb \
  --variants scratch/build1/variants_hg19_eaf.tsv \
  --traits   scratch/build1/traits.tsv \
  --variants-build hg19 \
  --ld-dir /local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38 \
  --ld-ancestry EUR \
  --workers 128 \
  --ld-vcf-threads 8 \
  --chunk-v 512 --chunk-t 512 \
  --overwrite > scratch/build3/main_build.log 2>&1 &
echo "PID: $!"
```

If interrupted, resume from checkpoint (no work is lost beyond the current
512-trait batch):

```bash
cd /home/gh13047/repo/pleiodb

pleiodb build \
  /local-scratch/data/pleiodb/main.pleiodb \
  --variants scratch/build1/variants_hg19_eaf.tsv \
  --traits   scratch/build1/traits.tsv \
  --variants-build hg19 \
  --ld-dir /local-scratch/projects/genotype-phenotype-map/data/ld_reference_panel_hg38 \
  --ld-ancestry EUR \
  --workers 128 \
  --ld-vcf-threads 8 \
  --chunk-v 512 --chunk-t 512 \
  --resume > scratch/build3/main_build.log 2>&1 &
echo "PID: $!"
```

Check progress:
```bash
tail -5 scratch/build3/main_build.log
python3 -c "import json; d=json.load(open('/local-scratch/data/pleiodb/main.pleiodb/build_checkpoint.json')); print(d['phase'], d['ingestion_t_batches_done'], '/', 9, 'ingestion |', d['imputation_t_batches_done'], '/', 9, 'imputation')"
```

## Rho Matrix (after build completes)

```bash
pleiodb rho \
  /local-scratch/data/pleiodb/main.pleiodb \
  --workers 32 2>&1 | tee scratch/build3/main_rho.log
```

## Verification

```bash
pleiodb info /local-scratch/data/pleiodb/main.pleiodb

# Imputation coverage
python3 -c "
import sys; sys.path.insert(0, 'src')
import numpy as np
from pleiodb.db import GWASDatabase
from pleiodb.quantize import decode_z
db = GWASDatabase.open('/local-scratch/data/pleiodb/main.pleiodb')
V, T = db.V, db.T
z = decode_z(db._zscore.get_block(0, V, 0, T))
n = V * T
fin = int(np.isfinite(z).sum())
print(f'V={V} T={T} total={n:,}')
print(f'Finite:  {fin:,}  ({100*fin/n:.2f}%)')
print(f'Missing: {n-fin:,}  ({100*(n-fin)/n:.2f}%)')
"
```

## Files in this directory

| File | Description |
|------|-------------|
| `README.md` | This file |
| `main_build.log` | Build log (created when build runs) |
| `main_rho.log` | Rho matrix build log (created after build) |
