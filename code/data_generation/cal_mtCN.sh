
#!/bin/bash

module load bcftools/1.12

bcftools merge ~/mitonuclear2/data/dx_downloads/vcf/dump/*.vcf.gz \
-Oz -o ~/mitonuclear2/data/dx_downloads/vcf/pmbb_1_35000k.vcf.gz

find . -regextype posix-extended -regex ".*/pmbb_[3][5][0-9]{3}_[0-9]{5}.MT.vcf.gz" \
> ../mergelist.txt

cat ../mergelist.txt | while read filename; do \
bcftools index ${filename}; done

bcftools merge ~/mitonuclear2/data/dx_downloads/vcf/dump/pmbb_1_35000.vcf.gz \
~/mitonuclear2/data/dx_downloads/vcf/dump/pmbb_1_40000.vcf.gz

find . -regextype posix-extended -regex ".*/pmbb_[4][0-9]{4}_[0-9]{5}.MT.vcf.gz" \
> ../mergelist.txt

#add pmbb_1_40000.MT.vcf.gz to list

cat ../mergelist.txt | while read filename; do \
bcftools index ${filename}; done

bcftools merge -l ../mergelist.txt \
-Oz -o ~/mitonuclear2/data/dx_downloads/vcf/dump/pmbb_1_45000.vcf.gz

bcftools query -l ~/mitonuclear2/data/dx_downloads/vcf/dump/pmbb_1_45000.vcf.gz \
> ~/mitonuclear2/data/dx_downloads/vcf/pmbb_1_45000.samples

bcftools query -f '%POS[\t%DP]\n' pmbb_1_45000.vcf.gz > pmbb_1_45000.dp
gzip pmbb_1_45000.dp
