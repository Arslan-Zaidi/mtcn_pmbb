#!/bin/bash

module load gcta/1.93.2b
module load plink/2.0-20210505

#AST
mkdir -p ~/mitonuclear2/data/gcta/AFR/mtcn.experiments
gcta64 \
--grm ~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR \
--thread-num 11 \
--reml \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--qcovar ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/pmbb_AFR.pcs.ast_09022022.qcovar \
--out ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.ast

#DARC recessive
gcta64 \
--grm ~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR \
--thread-num 11 \
--reml \
--qcovar ~/mitonuclear2/data/gcta/AFR/pmbb_AFR.pcs.darc.a_09022022.qcovar \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--out ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.darc.a

gcta64 \
--grm ~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR \
--thread-num 11 \
--reml \
--qcovar ~/mitonuclear2/data/gcta/AFR/pmbb_AFR.pcs.darc.ad_09022022.qcovar \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--out  ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.darc.ad

#### split chromosome 1 into 4 parts
gcta64 \
--mgrm ~/mitonuclear2/data/gcta/AFR/mgrm.bychrom.chr1split.txt \
--thread-num 11 \
--reml \
--qcovar ~/mitonuclear2/data/gcta/AFR/pmbb_AFR.pcs.darc.ad_09022022.qcovar \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--out ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.darc.ad.chrom1

#### APOL1
gcta64 \
--grm ~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR \
--thread-num 11 \
--reml \
--qcovar ~/mitonuclear2/data/gcta/AFR/pmbb_AFR.pcs.apol1_09062022.qcovar \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--out  ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.apol1


##### admixture mapping hits

gcta64 \
--grm ~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR \
--thread-num 11 \
--reml \
--qcovar ~/mitonuclear2/data/gcta/AFR/pmbb_AFR.pcs.lanchits_09062022.qcovar \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--out  ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.lanchits


#### blood prs

gcta64 \
--grm ~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR \
--thread-num 11 \
--reml \
--qcovar ~/mitonuclear2/data/gcta/AFR/pmbb_AFR.pcs.blprs_09082022.qcovar \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--out  ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.blprs

gcta64 \
--grm ~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR \
--thread-num 11 \
--reml \
--qcovar ~/mitonuclear2/data/gcta/AFR/pmbb_AFR.pcs.blprs_09082022.extended.qcovar \
--pheno ~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno \
--out  ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.blprs.ext


#concatenate all together
alt_models=("apol1" "ast" "blprs.ext" "blprs" "darc.a" "darc.ad" "lanchits")
for i in ${alt_models[@]}; do \
grep "V(G)/Vp" ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.${i}.hsq | \
awk -v OFS="\t" -v alt=$i '{print alt,$0}' ; done \
> ~/mitonuclear2/data/gcta/AFR/mtcn.experiments/PMBB_reml_AFR.mtcn.logmtresid.all.hsq
