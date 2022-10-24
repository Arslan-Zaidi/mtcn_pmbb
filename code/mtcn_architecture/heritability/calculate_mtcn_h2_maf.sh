#!/bin/bash

module load gcta/1.93.2b
#module load plink/2.0-20210505

pop=${1}
covar=${2}

if [ "$covar" == "nocovar" ]; then
#w/o covariates
gcta64 \
--thread-num 11 \
--mgrm ~/mitonuclear2/data/gcta/${pop}/mgrm.bymaf.txt \
--reml \
--pheno ~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_07122022.txt \
--mpheno 2 \
--out ~/mitonuclear2/data/gcta/${pop}/maf/PMBB_reml_${pop}.logmtcn.nocovar
else
  if [ "$covar" == "blood" ]; then
    #covariates
    gcta64 \
    --mgrm ~/mitonuclear2/data/gcta/${pop}/mgrm.bymaf.txt \
    --thread-num 11 \
    --reml \
    --covar ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.sex \
    --pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_residuals_07192022.pheno \
    --out ~/mitonuclear2/data/gcta/${pop}/maf/PMBB_reml_${pop}.logmtcn.blresid
  else
    if [ "$covar" == "DARC" ]; then
      gcta64 \
      --mgrm ~/mitonuclear2/data/gcta/${pop}/mgrm.bymaf.txt \
      --thread-num 11 \
      --reml \
      --covar ~/mitonuclear2/data/genotype/DARC/pmbb_rs281477.${pop}.raw2 \
      --pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_residuals_07192022.pheno \
      --out ~/mitonuclear2/data/gcta/${pop}/maf/PMBB_reml_${pop}.logmtcn.blresid.darc
    fi
  fi
fi

