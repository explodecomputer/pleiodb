#!/bin/bash

conda activate bcftools

# I have a list of variants in format 1:1231231_G_T
# Convert it to format 1\t1231231\t1231231

igddir="/local-scratch/data/opengwas/igd"
mkdir -p tests/test_data/vcf
trait_list=$(cut -f 1 tests/test_data/traits.tsv)

parse_variants() {
    input=$1
    output=$2
    tmp1=$(mktemp)
    cut -d ":" -f 1 $input > $tmp1
    cut -d ":" -f 2 $input | cut -d "_" -f 1 > $tmp1.pos
    paste -d "\t" $tmp1 $tmp1.pos > $output
}

parse_variants tests/test_data/input_variants_hg19.txt tests/test_data/input_variants_hg19_parsed.txt

for trait in $trait_list; do
    echo "Processing trait: $trait"
    bcftools view -R tests/test_data/input_variants_hg19_parsed.txt -Oz $igddir/$trait/$trait.vcf.gz > tests/test_data/vcf/$trait.vcf.gz
    bcftools index --force tests/test_data/vcf/$trait.vcf.gz
    bcftools index -n tests/test_data/vcf/$trait.vcf.gz
done

bcftools query -f '%CHROM:%POS[_]%REF[_]%ALT\t%AF\n' tests/test_data/vcf/ieu-a-7.vcf.gz | shuf | head -n 500 | sort > tests/test_data/variants_hg19.tsv


