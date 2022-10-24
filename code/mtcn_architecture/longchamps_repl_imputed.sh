#!/bin/bash

module load plink/2.0-20210505

pop=${1}
chrom=${2}

mkdir -p ~/mitonuclear2/data/mtgwas/mtCN/imputed/lchamps_rep/${pop}

#filtered log mtCN

#with blood traits as covariates

plink2 \
--pfile ~/pmbb/Imputed/pgen/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_chr${chrom} \
--keep ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_07152022.covar \
--glm hide-covar \
--extract ~/mitonuclear2/data/prev_associations/lchamps_pmbbid.txt \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--covar ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_07152022.covar \
--covar-col-nums 4-23 \
--out ~/mitonuclear2/data/mtgwas/mtCN/imputed/lchamps_rep/${pop}/residuals/pmbb_lchamps_logmtcn.resid.imputed.${pop}.chr${chrom}
