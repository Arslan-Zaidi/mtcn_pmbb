#!/bin/bash

module load gcta/1.93.2b
#module load plink/2.0-20210505

pop=${1}

gcta64 \
--mgrm ~/mitonuclear2/data/gcta/${pop}/mgrm.bychrom.txt \
--thread-num 11 \
--reml \
--covar ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.sex \
--qcovar ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.noblood.qcovar \
--pheno ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.blood.qcovar \
--mpheno 27 \
--out ~/mitonuclear2/data/gcta/${pop}/chrom/PMBB_reml_${pop}.neutrophil.covar

gcta64 \
--mgrm ~/mitonuclear2/data/gcta/${pop}/mgrm.bychrom.txt \
--thread-num 11 \
--reml \
--covar ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.sex \
--qcovar ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.noblood.rs2814778.qcovar \
--pheno ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.blood.qcovar \
--mpheno 27 \
--out ~/mitonuclear2/data/gcta/${pop}/chrom/PMBB_reml_${pop}.neutrophil.rs2814778.covar
