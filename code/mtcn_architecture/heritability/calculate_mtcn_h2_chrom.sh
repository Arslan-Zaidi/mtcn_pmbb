#!/bin/bash

module load gcta/1.93.2b
#module load plink/2.0-20210505

pop=${1}
covar=${2}

# #split covar and qcovar
# cut -f1-3 ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_03042022.covar \
# > ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.sex
#
# #split covar and qcovar
# cut -f1,2,4-25 ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_03042022.covar \
# > ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.noblood.qcovar

if [ "$covar" == "nocovar" ]; then
#w/o covariates
gcta64 \
--thread-num 11 \
--mgrm ~/mitonuclear2/data/gcta/${pop}/mgrm.bychrom.txt \
--reml \
--pheno ~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_07122022.txt \
--mpheno 2 \
--out ~/mitonuclear2/data/gcta/${pop}/chrom/PMBB_reml_${pop}.logmtcn.nocovar
else
  if [ "$covar" == "pcs" ]; then
    #covariates
    gcta64 \
    --mgrm ~/mitonuclear2/data/gcta/${pop}/mgrm.bychrom.txt \
    --thread-num 11 \
    --reml \
    --pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
    --qcovar ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop}.maf0.01.pca.header.eigenvec \
    --out ~/mitonuclear2/data/gcta/${pop}/chrom/PMBB_reml_${pop}.logmtcn.blresid
  else
    if [ "$covar" == "DARC" ]; then
      gcta64 \
      --mgrm ~/mitonuclear2/data/gcta/${pop}/mgrm.bychrom.txt \
      --thread-num 11 \
      --reml \
      --qcovar ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.pcs.rs2814778.qcovar \
      --pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
      --out ~/mitonuclear2/data/gcta/${pop}/chrom/PMBB_reml_${pop}.logmtcn.blresid.darc
    fi
  fi
fi
