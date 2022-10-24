#!/bin/bash

#partitioning heritability into mito genes and non-mito genes

pop=${1}

module load plink/2.0-20210505
module load gcta/1.93.2b

######## make global grm ################

#filter on maf and prune for LD

plink2 \
--bfile ~/pmbb/Genotype/PMBB-Release-2020-2.0_genetic_genotype \
--chr 1-22 \
--keep ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_03042022.covar \
--make-bed \
--out ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.autosome.${pop}

plink2 \
--bfile ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.autosome.${pop} \
--maf 0.01 \
--snps-only just-acgt \
--indep-pairwise 100 10 0.1 \
--out ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.${pop}

plink2 \
--bfile ~/pmbb/Genotype/PMBB-Release-2020-2.0_genetic_genotype \
--keep ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_03042022.covar \
--extract ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.${pop}.prune.in \
--make-bed --out ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.${pop}

gcta64 \
--bfile ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.${pop} \
--make-grm --out ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop}

######## remove rarer variants

plink2 \
--bfile ~/pmbb/Genotype/PMBB-Release-2020-2.0_genetic_genotype \
--chr 1-22 \
--keep ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_03042022.covar \
--make-bed \
--out ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.autosome.${pop}

plink2 \
--bfile ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.autosome.${pop} \
--maf 0.05 \
--snps-only just-acgt \
--indep-pairwise 100 10 0.1 \
--out ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.05.acgt.ldpruned.${pop}

plink2 \
--bfile ~/pmbb/Genotype/PMBB-Release-2020-2.0_genetic_genotype \
--keep ~/mitonuclear2/data/phenotype/PMBB_plink4mat_${pop}.blood_03042022.covar \
--extract ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.05.acgt.ldpruned.${pop}.prune.in \
--make-bed --out ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.05.acgt.ldpruned.${pop}

gcta64 \
--bfile ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.05.acgt.ldpruned.${pop} \
--make-grm --out ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop}.maf0.05


#####  partition grm into mito and non-mito genes #####
rm ~/mitonuclear2/data/gcta/${pop}/mgrm1.txt
mkdir -p ~/mitonuclear2/data/gcta/${pop}/mtcarta/

plink2 --bfile ~/pmbb/Genotype/PMBB-Release-2020-2.0_genetic_genotype \
--extract bed0 ~/mitonuclear2/data/mitocarta3_0based_hg38.txt \
--maf 0.01 \
--snps-only just-acgt \
--keep ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.${pop}.fam \
--make-bed \
--out ~/mitonuclear2/data/gcta/${pop}/mtcarta/PMBB-Release-2020-2.0_genetic_maf0.01.acgt.${pop}_mtcarta

plink2 --bfile ~/mitonuclear2/data/gcta/${pop}/mtcarta/PMBB-Release-2020-2.0_genetic_maf0.01.acgt.${pop}_mtcarta \
--indep-pairwise 100 10 0.1 \
--out ~/mitonuclear2/data/gcta/${pop}/mtcarta/PMBB-Release-2020-2.0_genetic_maf0.01.acgt.ldpruned.${pop}_mtcarta

plink2 --bfile ~/mitonuclear2/data/gcta/${pop}/mtcarta/PMBB-Release-2020-2.0_genetic_maf0.01.acgt.${pop}_mtcarta \
--extract ~/mitonuclear2/data/gcta/${pop}/mtcarta/PMBB-Release-2020-2.0_genetic_maf0.01.acgt.ldpruned.${pop}_mtcarta.prune.in \
--make-bed \
--out ~/mitonuclear2/data/gcta/${pop}/mtcarta/PMBB-Release-2020-2.0_genetic_maf0.01.acgt.ldpruned.${pop}_mtcarta

gcta64 \
--bfile  ~/mitonuclear2/data/gcta/${pop}/mtcarta/PMBB-Release-2020-2.0_genetic_maf0.01.acgt.ldpruned.${pop}_mtcarta \
--make-grm --out ~/mitonuclear2/data/gcta/${pop}/mtcarta/PMBB-Release-2020-2.0_genetic.${pop}_mtcarta

touch ~/mitonuclear2/data/gcta/${pop}/mgrm1.txt
echo ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.${pop} >> \
~/mitonuclear2/data/gcta/${pop}/mgrm1.txt
echo ~/mitonuclear2/data/gcta/${pop}/mtcarta/PMBB-Release-2020-2.0_genetic.${pop}_mtcarta >> \
~/mitonuclear2/data/gcta/${pop}/mgrm1.txt


###### partition grm into chromosomes #########
rm ~/mitonuclear2/data/gcta/${pop}/mgrm.bychrom.txt
touch ~/mitonuclear2/data/gcta/${pop}/mgrm.bychrom.txt
mkdir -p ~/mitonuclear2/data/gcta/${pop}/chrom/

for chrom in {1..22}; do \
plink2 \
--bfile ~/mitonuclear2/data/gcta/${pop}/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.${pop} \
--chr ${chrom} \
--make-bed --out ~/mitonuclear2/data/gcta/${pop}/chrom/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.${pop}.${chrom}; \
done

for chrom in {1..22}; do \
gcta64 \
--bfile ~/mitonuclear2/data/gcta/${pop}/chrom/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.${pop}.${chrom} \
--make-grm --out ~/mitonuclear2/data/gcta/${pop}/chrom/PMBB-Release-2020-2.0_genetic_genotype.${pop}.${chrom}; \
done

for chrom in {1..22}; do \
echo ~/mitonuclear2/data/gcta/${pop}/chrom/PMBB-Release-2020-2.0_genetic_genotype.${pop}.${chrom} >> \
~/mitonuclear2/data/gcta/${pop}/mgrm.bychrom.txt; \
done


######## partition chromosome 1 into four
tail -n+2 ~/mitonuclear2/data/gcta/AFR/mgrm.bychrom.txt > \
~/mitonuclear2/data/gcta/AFR/mgrm.bychrom.chr1split.txt


47632625
111365828
203997762
248914573

count=1
cat ~/mitonuclear2/data/gcta/AFR/chr1.quart.bed | \
while read chrom start end; do \
plink2 --bfile ~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.AFR \
--chr ${chrom} --from-bp ${start} --to-bp ${end} \
--make-bed --out ~/mitonuclear2/data/gcta/AFR/chrom/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.AFR.${chrom}_${count}; \
count=$(( count+1 )); \
done

for i in {1..4}; do \
gcta64 \
--bfile ~/mitonuclear2/data/gcta/AFR/chrom/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.AFR.chr1_${i} \
--make-grm --out ~/mitonuclear2/data/gcta/AFR/chrom/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.AFR.chr1_${i}; \
done

for i in {1..4}; do
echo ~/mitonuclear2/data/gcta/AFR/chrom/PMBB-Release-2020-2.0_genetic_genotype.maf0.01.acgt.ldpruned.AFR.chr1_${i}; \
done >> ~/mitonuclear2/data/gcta/AFR/mgrm.bychrom.chr1split.txt
