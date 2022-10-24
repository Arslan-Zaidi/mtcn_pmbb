#!/bin/bash

#format: PMBB_plink4mat_${pop}.${covar}.qcovar

# #split covar and qcovar
#covar
# cut -f1-3 ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_07152022.covar \
# > ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.sex

#qcovar
# no blood traits
# cut -f1,2,4-23 ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_07152022.covar \
# > ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.noblood.qcovar

# with blood traits
# cut -f1,2,4-25,29-34 ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_07152022.covar \
# > ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.blood.qcovar

cut -f1,2,7 ~/mitonuclear2/data/genotype/DARC/pmbb_rs281477.${pop}.raw > \
~/mitonuclear2/data/genotype/DARC/pmbb_rs281477.${pop}.raw2

# join -j1 \
# <(sort -k1 ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.blood.qcovar) \
# <(sort -k1 ~/mitonuclear2/data/genotype/DARC/pmbb_rs281477.${pop}.raw2) \
# > ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.blood.rs2814778.qcovar

# no blood traits but DARC locus (for neutrophil h2)
join -j1 \
<(sort -k1 ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop}.maf0.01.pca.header.eigenvec) \
<(cut -f 2,3 ~/mitonuclear2/data/genotype/DARC/pmbb_rs281477.${pop}.raw2 | sort -k1) \
> ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.pcs.rs2814778.qcovar

#only sex and DARC (for resid logmtcn GWAS)
# join -j1 \
# <(sort ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.sex) \
# <(sort -k1 ~/mitonuclear2/data/genotype/DARC/pmbb_rs281477.${pop}.raw2) \
# > ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.sex.rs2814778.qcovar

# with blood traits (maf 0.05)
cat \
<(head -n1 ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_07152022.covar | cut -f1,4-23) \
<(cat ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop}.maf0.05.pca.eigenvec) \
> ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop}.maf0.05.pca.header.eigenvec


join -j1 \
<(cut -f1,2,24,25,29-35 ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_07152022.covar | sort -k1) \
<(cut -d ' ' -f1,3-22 ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop}.maf0.05.pca.header.eigenvec | sort -k1) \
> ~/mitonuclear2/data/gcta/${pop}/PMBB_plink4mat_${pop}.blood.maf0.05.qcovar
