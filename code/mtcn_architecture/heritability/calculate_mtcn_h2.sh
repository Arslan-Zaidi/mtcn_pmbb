#!/bin/bash

module load gcta/1.93.2b
module load plink/2.0-20210505

pop=${1}


gcta64 \
--grm ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop} \
--thread-num 11 \
--reml \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--qcovar ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop}.maf0.01.pca.eigenvec \
--out ~/mitonuclear2/data/gcta/${pop}/PMBB_reml_${pop}.mtcn.logmtresid
