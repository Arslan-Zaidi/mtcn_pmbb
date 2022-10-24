#!/bin/bash

module load bcftools/1.12
module load plink/2.0-20210505
module load jre/1.8.0_211
module load python/2.7.12

source ~/.bash_profile

shopt -s expand_aliases
source ~/.bash_aliases

chrom=${1}

#mkdir -p ~/mitonuclear2/data/genotype/imputed/genotyped_snps/vcf
mkdir -p ~/mitonuclear2/data/rfmix/1kg_vcf
mkdir -p ~/mitonuclear2/data/rfmix/snplist
mkdir -p ~/mitonuclear2/data/rfmix/phased_vcf
mkdir -p ~/mitonuclear2/data/rfmix/merged_vcf
mkdir -p ~/mitonuclear2/data/rfmix/classes
mkdir -p ~/mitonuclear2/data/rfmix/output

#create a file with AFR cohort - needs to be fone only once so comment out in script
# tail -n+2 ~/mitonuclear2/data/phenotype/PMBB_plink4mat_AFR.blood_03042022.covar | \
# cut -f1 > ~/mitonuclear2/data/rfmix/pmbb.afr.samples

echo "Selecting variants on chromosome and AFR individuals in PMBB"

#splt vcf into chromosomes
bcftools view \
-r ${chrom} \
-v snps \
-m2 -M2 \
--force-samples \
-S ~/mitonuclear2/data/rfmix/pmbb.afr.samples \
-Oz -o ~/mitonuclear2/data/rfmix/vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.vcf.gz \
~/pmbb/Genotype/PMBB-Release-2020-2.0_genetic_genotype.vcf.gz

echo "removing duplicated variants from PMBB data"
bcftools index -f ~/mitonuclear2/data/rfmix/vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.vcf.gz

#vcf has 8955 snps

# bcftools view -H ~/mitonuclear2/data/rfmix/vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.vcf.gz | \
# cut -f3,1,2,4,5 | grep -v "GSA" > ~/mitonuclear2/data/rfmix/snplist/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.snplist

bcftools norm \
-d all \
-Oz -o ~/mitonuclear2/data/rfmix/vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.rmdup.vcf.gz \
~/mitonuclear2/data/rfmix/vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.vcf.gz

#vcf has 8938 snps now after removing duplicates

echo "subsetting 1kg data for CEU and YRI & making overlap"
#create list of variants in the phased data (to filter 1kg data on)
bcftools view -H \
~/mitonuclear2/data/rfmix/vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.rmdup.vcf.gz |
cut -f 1,2 > ~/mitonuclear2/data/rfmix/snplist/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.snplist

#extract these variants and CEU and YRI samples from 1kg data
bcftools view \
--force-samples \
-T ~/mitonuclear2/data/rfmix/snplist/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.snplist \
-v snps \
-S ~/mitonuclear2/data/rfmix/1kg_ceu.yri.list \
-Oz -o ~/mitonuclear2/data/rfmix/1kg_vcf/1kg.ceuryi.${chrom}.phased.vcf.gz \
/project/mathilab/data/1000g/phase3_vcfs/ALL.chr${chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

bcftools index -f ~/mitonuclear2/data/rfmix/1kg_vcf/1kg.ceuryi.${chrom}.phased.vcf.gz

#snps in 1kg: 8544

bcftools view -H \
~/mitonuclear2/data/rfmix/1kg_vcf/1kg.ceuryi.${chrom}.phased.vcf.gz | \
cut -f1,2 > ~/mitonuclear2/data/rfmix/snplist/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.snplist2

bcftools view \
-T ~/mitonuclear2/data/rfmix/snplist/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.snplist2 \
-Oz -o ~/mitonuclear2/data/rfmix/vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.rmdup.overlap.vcf.gz \
~/mitonuclear2/data/rfmix/vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.rmdup.vcf.gz

#snps in pmbb: 8544

#good to merge now

#create a list of snps to filter the pgen file for (archived)
# awk -v chr=${chrom} -v OFS="\t" \
# '$1==chr && NR==FNR {a[$4]; next} $2 in a {print $1,$2,$2,$3}' \
# ~/pmbb/Genotype/PMBB-Release-2020-2.0_genetic_genotype.bim \
# ~/pmbb/Imputed/pgen/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_chr${chrom}.pvar \
# > ~/mitonuclear2/data/genotype/imputed/genotyped_snps/PMBB-Release-2020-2.0_genetic_genotype.${chrom}.snplist

echo "running beagle -phasing and imputing PMBB data"
#run beagle
beagle gt=~/mitonuclear2/data/rfmix/vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.rmdup.overlap.vcf.gz \
out=~/mitonuclear2/data/rfmix/phased_vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.phased \
map=/project/mathilab/data/maps/beagle_maps/plink.chr${chrom}.GRCh38.map

#no. of snps: 8544 (good!)

echo "merging PMBB and 1kg VCFs"

bcftools index -f ~/mitonuclear2/data/rfmix/1kg_vcf/1kg.ceuryi.${chrom}.phased.vcf.gz
bcftools index -f ~/mitonuclear2/data/rfmix/phased_vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.phased.vcf.gz

#merge with 1000 genomes
bcftools merge \
-m both \
-Oz -o ~/mitonuclear2/data/rfmix/merged_vcf/pmbb.afr.1kg.ceuryi.${chrom}.phased.vcf.gz \
~/mitonuclear2/data/rfmix/phased_vcf/PMBB-Release-2020-2.0_genetic_genotype.AFR.${chrom}.phased.vcf.gz \
~/mitonuclear2/data/rfmix/1kg_vcf/1kg.ceuryi.${chrom}.phased.vcf.gz

#remove multi-allelic records
bcftools view \
-v snps \
-m2 -M 2 \
-o ~/mitonuclear2/data/rfmix/merged_vcf/pmbb.afr.1kg.ceuryi.${chrom}.phased.biallelic.vcf.gz -Oz \
~/mitonuclear2/data/rfmix/merged_vcf/pmbb.afr.1kg.ceuryi.${chrom}.phased.vcf.gz

echo "creating RFMIX alleles file"
#Generating input files for rfmix
#1. alleles file
#the rfmix format has snps in rows and haplotypes in columns without spaces
bcftools query -f '[%GT]\n' \
~/mitonuclear2/data/rfmix/merged_vcf/pmbb.afr.1kg.ceuryi.${chrom}.phased.biallelic.vcf.gz | \
tr -d '|' \
> ~/mitonuclear2/data/rfmix/alleles/pmbb.afr.1kg.ceuryi.${chrom}.phased.alleles

echo "creating SNP files"
#2. snp_locations file
#genetic coordinates in cM in the same order as rows in the alleles file

bcftools view -H \
~/mitonuclear2/data/rfmix/merged_vcf/pmbb.afr.1kg.ceuryi.${chrom}.phased.biallelic.vcf.gz | \
cut -f1,2 > ~/mitonuclear2/data/rfmix/snplist/pmbb.afr.1kg.ceuryi.${chrom}.phased.snplist

Rscript ~/mitonuclear2/code/cal_gmap.R \
-g /project/mathilab/data/maps/beagle_maps/plink.chr${chrom}.GRCh38.map \
-b ~/mitonuclear2/data/rfmix/snplist/pmbb.afr.1kg.ceuryi.${chrom}.phased.snplist \
-o ~/mitonuclear2/data/rfmix/snplist/pmbb.afr.1kg.ceuryi.${chrom}.phased.snplist

cut -f3 \
~/mitonuclear2/data/rfmix/snplist/pmbb.afr.1kg.ceuryi.${chrom}.phased.snplist \
> ~/mitonuclear2/data/rfmix/snp_locations/pmbb.afr.1kg.ceuryi.${chrom}.phased.snp_locations

echo "creating samples file"
bcftools query -l \
~/mitonuclear2/data/rfmix/merged_vcf/pmbb.afr.1kg.ceuryi.${chrom}.phased.vcf.gz \
> ~/mitonuclear2/data/rfmix/classes/pmbb.afr.1kg.ceuryi.${chrom}.phased.ids

Rscript ~/mitonuclear2/code/generate_rfmix_classes.R \
-i ~/mitonuclear2/data/rfmix/classes/pmbb.afr.1kg.ceuryi.${chrom}.phased.ids \
-p ~/mitonuclear2/data/rfmix/1kg_ceu.yri.pop \
-o ~/mitonuclear2/data/rfmix/classes/pmbb.afr.1kg.ceuryi.${chrom}.phased.classes


echo "Running rfmix"
python RunRFMix.py \
--forward-backward \
PopPhased \
-o ~/mitonuclear2/data/rfmix/output/pmbb.afr.1kg.ceuryi.${chrom} \
~/mitonuclear2/data/rfmix/alleles/pmbb.afr.1kg.ceuryi.${chrom}.phased.alleles \
~/mitonuclear2/data/rfmix/classes/pmbb.afr.1kg.ceuryi.${chrom}.phased.classes \
~/mitonuclear2/data/rfmix/snp_locations/pmbb.afr.1kg.ceuryi.${chrom}.phased.snp_locations
