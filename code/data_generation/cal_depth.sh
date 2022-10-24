#!/bin/bash

cat ~/mitonuclear2/data/genotype/PMBB_Release-2020-2.0-MT.sites.subsampled.from.exome.GL.txt | \
while read chrom pos; do \
bcftools query -r ${chrom}:${pos} -f '%CHROM\t%POS\t[\t%DP]\n' ~/pmbb/Exome/pVCF/GL_by_chrom/PMBB-Release-2020-2.0_genetic_exome_chr${chrom}_GL.vcf.gz; \
done > ~/mitonuclear2/data/genotype/PMBB_Release-2020-2.0-MT.sites.subsampled.from.exome.GL.DP


bcftools query -l ~/pmbb/Exome/pVCF/GL_by_chrom/PMBB-Release-2020-2.0_genetic_exome_chr22_GL.vcf.gz \
> ~/mitonuclear2/data/genotype/PMBB_Release-2020-2.0-MT.sites.subsampled.from.exome.GL.samples
