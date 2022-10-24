#!/bin/bash

module load plink/2.0-20210505

pop=${1}
chrom=${2}

mkdir -p ~/mitonuclear2/data/mtgwas/mtCN/imputed/${pop}

join -j1 \
<(sort -k1 ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.blood.rs2814778.qcovar) \
<(cut -f1,3 ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.sex | sort -k1) \
> ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood.rs2814778.covar

#filtered log mtCN

#with blood traits as covariates

plink2 \
--pfile ~/pmbb/Imputed/pgen/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_chr${chrom} \
--glm hide-covar \
--keep ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood.rs2814778.covar \
--maf 0.01 \
--pheno ~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_07122022.txt \
--pheno-col-nums 4 \
--covar ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood.rs2814778.covar \
--out ~/mitonuclear2/data/mtgwas/mtCN/imputed/${pop}/pmbb_mtgwas_logmtcn.filtered.imputed.${pop}.chr${chrom}.bl.rs2814778 \
--covar-col-nums 3-33
