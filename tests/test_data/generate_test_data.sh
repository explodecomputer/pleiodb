#!/bin/bash
# Regenerate test data from the OpenGWAS IGD VCF archive.
#
# Prerequisites: activate the pleiodb conda environment first.
#   mamba activate pleiodb
#
# Usage (from repo root):
#   bash tests/test_data/generate_test_data.sh
#
# Outputs:
#   tests/test_data/vcf/            — one bgzipped VCF + CSI index per trait
#   tests/test_data/variants_hg19.tsv  — 500 ALIDs + EAF, hg19
#   tests/test_data/variants_hg38.tsv  — same 500 variants lifted to hg38

set -euo pipefail

igddir="/local-scratch/data/opengwas/igd"
mkdir -p tests/test_data/vcf

trait_list=$(cut -f 1 tests/test_data/traits.tsv)

# ---------------------------------------------------------------------------
# Parse input_variants_hg19.txt (CHROM:POS_REF_ALT) to CHROM\tPOS for bcftools
# ---------------------------------------------------------------------------
parse_variants() {
    input=$1
    output=$2
    tmp1=$(mktemp)
    cut -d ":" -f 1 "$input" > "$tmp1"
    cut -d ":" -f 2 "$input" | cut -d "_" -f 1 > "$tmp1.pos"
    paste -d $'\t' "$tmp1" "$tmp1.pos" > "$output"
    rm -f "$tmp1" "$tmp1.pos"
}

parse_variants tests/test_data/input_variants_hg19.txt \
               tests/test_data/input_variants_hg19_parsed.txt

# ---------------------------------------------------------------------------
# Extract per-trait VCF subsets
# ---------------------------------------------------------------------------
for trait in $trait_list; do
    vcf_path="$igddir/$trait/$trait.vcf.gz"
    out_path="tests/test_data/vcf/$trait.vcf.gz"
    echo "Processing trait: $trait"
    bcftools view -R tests/test_data/input_variants_hg19_parsed.txt \
                  -Oz "$vcf_path" > "$out_path"
    bcftools index --force "$out_path"
done

# ---------------------------------------------------------------------------
# Generate variants_hg19.tsv: 500 ALIDs + EAF in canonical allele order
#
# Canonical ALID: CHROM:POS_A1_A2  where A1 <= A2 alphabetically.
# EAF = frequency of A2 (the effect / alphabetically-second allele).
#
# For each record:
#   If REF <= ALT alphabetically: ALID = CHROM:POS_REF_ALT,  EAF = AF
#   If REF  > ALT alphabetically: ALID = CHROM:POS_ALT_REF,  EAF = 1 - AF
# ---------------------------------------------------------------------------
echo "Generating variants_hg19.tsv ..."

bcftools query -f '%CHROM:%POS\t%REF\t%ALT\t%AF\n' \
    tests/test_data/vcf/ieu-a-7.vcf.gz \
    | shuf | head -n 500 | sort \
    | python3 - <<'PYEOF'
import sys, math

lines = sys.stdin.read().splitlines()
rows = []
for line in lines:
    if not line.strip():
        continue
    chrom_pos, ref, alt, af_str = line.split('\t')
    try:
        af = float(af_str)
    except ValueError:
        af = math.nan
    if ref <= alt:
        alid = f"{chrom_pos}_{ref}_{alt}"
        eaf = af
    else:
        alid = f"{chrom_pos}_{alt}_{ref}"
        eaf = 1.0 - af if not math.isnan(af) else math.nan
    rows.append((alid, eaf))

for alid, eaf in rows:
    if math.isnan(eaf):
        print(alid)
    else:
        print(f"{alid}\t{eaf:.6g}")
PYEOF
# The python heredoc above reads from the pipe; redirect properly:
bcftools query -f '%CHROM:%POS\t%REF\t%ALT\t%AF\n' \
    tests/test_data/vcf/ieu-a-7.vcf.gz \
    | shuf | head -n 500 | sort \
    | python3 -c "
import sys, math
for line in sys.stdin:
    line = line.rstrip()
    if not line:
        continue
    chrom_pos, ref, alt, af_str = line.split('\t')
    try:
        af = float(af_str)
    except ValueError:
        af = math.nan
    if ref <= alt:
        alid = f'{chrom_pos}_{ref}_{alt}'
        eaf = af
    else:
        alid = f'{chrom_pos}_{alt}_{ref}'
        eaf = 1.0 - af if not math.isnan(af) else math.nan
    if math.isnan(eaf):
        print(alid)
    else:
        print(f'{alid}\t{eaf:.6g}')
" > tests/test_data/variants_hg19.tsv

echo "  → $(wc -l < tests/test_data/variants_hg19.tsv) variants written to variants_hg19.tsv"

# ---------------------------------------------------------------------------
# Generate variants_hg38.tsv: lift hg19 ALIDs to hg38 coordinates
# Alleles and EAF are unchanged (only coordinates shift).
# ---------------------------------------------------------------------------
echo "Lifting variants_hg19.tsv → variants_hg38.tsv ..."

python3 -c "
import sys
from pyliftover import LiftOver

lo = LiftOver('hg19', 'hg38')
lifted = 0
failed = 0

with open('tests/test_data/variants_hg19.tsv') as fh:
    for line in fh:
        line = line.rstrip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        alid = parts[0]
        eaf_str = parts[1] if len(parts) > 1 else ''

        # Parse ALID
        colon = alid.index(':')
        chrom = alid[:colon]
        rest = alid[colon+1:]
        pos_str, a1, a2 = rest.split('_', 2)
        pos = int(pos_str)

        chrom_in = chrom if chrom.startswith('chr') else f'chr{chrom}'
        result = lo.convert_coordinate(chrom_in, pos - 1)
        if not result:
            failed += 1
            continue

        new_chrom = result[0][0].lstrip('chr')
        new_pos = int(result[0][1]) + 1
        new_alid = f'{new_chrom}:{new_pos}_{a1}_{a2}'

        if eaf_str:
            print(f'{new_alid}\t{eaf_str}')
        else:
            print(new_alid)
        lifted += 1

print(f'Lifted {lifted} variants; {failed} failed', file=sys.stderr)
" > tests/test_data/variants_hg38.tsv

echo "  → $(wc -l < tests/test_data/variants_hg38.tsv) variants written to variants_hg38.tsv"
echo "Done."
