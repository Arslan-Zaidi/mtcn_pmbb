#!/bin/bash

module load plink/2.0-20210505

pop=${1}
chrom=${2}

mkdir -p ~/mitonuclear2/data/mtgwas/mtCN/imputed/${pop}

# plink2 \
# --bfile ~/pmbb/Genotype/PMBB-Release-2020-2.0_genetic_genotype \
# --glm hide-covar \
# --chr 1-22,26 \
# --pheno ~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_02092022.txt \
# --covar ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.covar \
# --out ~/mitonuclear2/data/mtgwas/mtCN/${pop}/pmbb_mtgwas_mtCN.filtered.${pop}.nobl \
# --covar-col-nums 3-25
#
# #with blood traits as covariates
#
# plink2 \
# --bfile ~/pmbb/Genotype/PMBB-Release-2020-2.0_genetic_genotype \
# --glm hide-covar \
# --chr 1-22,26 \
# --pheno ~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_02092022.txt \
# --covar ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood.covar \
# --out ~/mitonuclear2/data/mtgwas/mtCN/${pop}/pmbb_mtgwas_mtCN.filtered.${pop}.bl \
# --covar-col-nums 3-25,29-34

#filtered INT

plink2 \
--pfile ~/pmbb/Imputed/pgen/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_chr${chrom} \
--glm hide-covar \
--maf 0.01 \
--snps-only \
--pheno ~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_02092022.txt \
--covar ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_03042022.covar \
--out ~/mitonuclear2/data/mtgwas/mtCN/imputed/${pop}/pmbb_mtgwas_mtCN.filtered.imputed.${pop}.chr${chrom}.nobl \
--covar-col-nums 3-25

#with blood traits as covariates

plink2 \
--pfile ~/pmbb/Imputed/pgen/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_chr${chrom} \
--glm hide-covar \
--maf 0.01 \
--snps-only \
--pheno ~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_02092022.txt \
--covar ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_03042022.covar \
--out ~/mitonuclear2/data/mtgwas/mtCN/imputed/${pop}/pmbb_mtgwas_mtCN.filtered.imputed.${pop}.chr${chrom}.bl \
--covar-col-nums 3-25,29-35
