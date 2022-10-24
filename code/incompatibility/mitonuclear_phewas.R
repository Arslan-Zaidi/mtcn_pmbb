#!/usr/bin/env Rscript

# author: AAZaidi

# suppressPackageStartupMessages(library('optparse')) option_list <- list(
# make_option(c('-g', '--glanc'), type = 'character', help = 'path to mean
# local ancestry in some regions', default = 'None'), make_option(c('-o',
# '--output'), type = 'character', help = 'path of output file'),
# make_option(c('-p', '--phecode1'), type = 'character', help = 'phecode for
# the phenotype to be tested, leave blank if testing against all phenotypes',
# default = 'None') ) args <- parse_args(OptionParser(option_list =
# option_list))

"%ni%" <- Negate("%in%")


suppressWarnings(suppressMessages({
  library(data.table)
  library(PheWAS)
  library(pbapply)
}))

# define root directory F = rprojroot::is_rstudio_project$make_fix_file()


# read ancestry information
ancestry = fread("~/mitonuclear2/data/genotype/pmbb.afr.weighted.04222022.glanc")

mtcn = fread("~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_10212022.pheno")



# if (args$glanc != 'None') { \t#read local ancestry from mitochondrial genes
# \tglanc = fread(args$glanc, header = FALSE) \tcolnames(glanc) = c('IID',
# 'glanc') \thaplo = merge(haplo, glanc, by = 'IID') \t#calculate correlation
# between global_ancestry and mean regoin_specific ancestry \tcor.glanc =
# with(haplo, cor(RFMIX, glanc)) \tanc_variable = 'glanc' }else{ \tanc_variable
# = 'RFMIX' }

phecodes = fread("~/pmbb/Phenotype/2.2/PMBB-Release-2020-2.2_phenotype_PheCode-matrix.txt",
  header = TRUE)
# phecodes.afr =
# fread('~/mitonuclear2/data/phenotype/PMBB_plink4mat_phecodes_AFR.pheno',
# header = TRUE) phecodes.afr = melt(phecodes.afr, id.vars = c('FID','IID'))
# phecodes.afr[value == -9, value := NA] phecodes.afr = dcast(phecodes.afr, FID
# + IID ~ variable)

bin.phenotypes = setdiff(colnames(phecodes), c("PMBB_ID"))

afr.covar = fread("~/mitonuclear2/data/phenotype/PMBB_plink4mat_AFR_07152022.covar")
afr.covar = afr.covar[IID %in% mtcn$IID,]
afr.covar[Haplo1 == "L", `:=`(regional_haplo, "African")]
afr.covar[Haplo1 %in% c("H", "I", "J", "K", "N", "R", "T", "U", "V", "W", "X"), `:=`(regional_haplo,
  "European")]
afr.covar = afr.covar[regional_haplo %in% c("African", "European")]
afr.covar$regional_haplo = factor(afr.covar$regional_haplo)
afr.covar$regional_haplo = relevel(afr.covar$regional_haplo, ref = "African")

afr.covar = melt(afr.covar, id.vars = c("FID", "IID", "Sex", "Haplogroup", "Haplo1",
  "Haplo2", "regional_haplo"))
afr.covar = afr.covar[-grep("PC", afr.covar$variable), ]
afr.covar[value == -9, `:=`(value, NA)]
afr.covar = dcast(afr.covar, FID + IID + Sex + Haplogroup + Haplo1 + Haplo2 + regional_haplo ~
  variable)
afr.covar = merge(afr.covar, ancestry[, .(IID, European)], by = "IID")

phecodes.afr = merge(afr.covar, phecodes, by.y = "PMBB_ID", by.x = "IID")




# function to run GWAS
logistic = function(pheno.code) {

  keepcols = c("IID", "regional_haplo", "Sex", "age_cent", "age_cent2", "European",
    pheno.code)

  dat = phecodes.afr[, keepcols, with = FALSE]
  colnames(dat)[7] = "nvisits"
  dat[nvisits >= 2, `:=`(phenotype, 1)]
  dat[nvisits == 0, `:=`(phenotype, 0)]
  # dat[case.control == 1, phenotype := 0] dat[case.control == 2, phenotype :=
  # 1]
  l1 = glm(data = dat, phenotype ~ regional_haplo + European + Sex + age_cent +
    age_cent2 + regional_haplo * European, family = binomial(link = "logit"),
    control = list(maxit = 1000))

  s1 = summary(l1)$coefficients
  ncases = length(which(dat$phenotype == 1))
  ncontrols = length(which(dat$phenotype == 0))

  return(data.table(phecode = rep(pheno.code, 7), predictor = c("Intercept", "European_Haplogroup",
    "European_ancestry", "Sex", "Age", "Age2", "Interaction"), beta = s1[, 1],
    zvalue = s1[, 3], pvalue = s1[, 4], convergence = l1$converged, ncases = ncases,
    ncontrols = ncontrols))

}

# carry out GWAS on all phenotypes
lresults1 = pblapply(bin.phenotypes, logistic)
dresults1 = dplyr::bind_rows(lresults1)

pheno.descriptions = addPhecodeInfo(unique(dresults1$phecode))
dresults1 = merge(dresults1, pheno.descriptions, by.x = "phecode", by.y = "phenotype")

# write this table
fwrite(dresults1, "~/mitonuclear2/data/mtgwas/mitonuc.incomp/pmbb_binary_phewas_mitonuc_10212022.txt",
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


# read median of the labs
mtcn1 = fread("~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_10212022.pheno")
mtcn1[, `:=`(variable, "rlrmtCN")]
mtcn1 = mtcn1[, .(IID, variable, rlrmtcn)]
colnames(mtcn1) = c("PMBB_ID", "Phenotype", "value")


quants = fread("~/mitonuclear2/data/phenotype/PMBB_2.2_quant_08122022.boxcox.pheno",
  header = TRUE)
colnames(quants)[3] = "value"

quants = rbind(mtcn1, quants)
quant.cols = unique(quants$Phenotype)

quants = merge(quants, afr.covar, by.x = "PMBB_ID", by.y = "IID")

# function to run GWAS
linear = function(phenotype) {

  # keepcols = c('PMBB_ID', 'regional_haplo', 'Sex', 'age_cent', 'age_cent2',
  # 'European', phenotype)

  dat = quants[Phenotype == phenotype, ]
  l1 = lm(data = dat, value ~ regional_haplo + European + Sex + age_cent + age_cent2 +
    regional_haplo * European)

  s1 = summary(l1)$coefficients
  nsamples = l1$df.residual + l1$rank

  return(data.table(phenotype = rep(phenotype, 7), predictor = c("Intercept", "European_Haplogroup",
    "European_ancestry", "Sex", "Age", "Age2", "Interaction"), beta = s1[, 1],
    zvalue = s1[, 3], pvalue = s1[, 4], nsamples = nsamples))

}

lm.afr.list = lapply(quant.cols, linear)
lm.afr.df = dplyr::bind_rows(lm.afr.list)

# write this table
fwrite(lm.afr.df, "~/mitonuclear2/data/mtgwas/mitonuc.incomp/pmbb_labs_phewas_mitonuc_bxcx_10212022.txt",
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


##### follow -up ####

## quantitative traits
qphenotypes = c("rlrmtCN", "CK", "Creatinine", "DIASTOLIC", "Neutrophil_median",
  "SYSTOLIC", "triglycerides")
qpheno.resids = lapply(qphenotypes, function(x) {
  dat = quants[Phenotype == x, .(PMBB_ID, Sex, age_cent, age_cent2, value, European,
    regional_haplo)]
  l1 = lm(data = dat, value ~ Sex + age_cent + age_cent2, na.action = "na.exclude")
  dat[, `:=`(residuals, resid(l1))]
  return(dat[, .(European, regional_haplo, residuals)])
})

names(qpheno.resids) = qphenotypes
qpheno.resids = dplyr::bind_rows(qpheno.resids, .id = "phenotype")


fwrite(qpheno.resids, "~/mitonuclear2/data/mtgwas/mitonuc.incomp/pmbb_qphenoresids_sex_age_age2_08142022.txt",
  sep = "\t")



## binary traits
bphenotypes = c("401", "401.1", "070.2")
bpheno.resids = lapply(bphenotypes, function(x) {
  keepcols = c("IID", "Sex", "age_cent", "age_cent2", x, "European", "regional_haplo")
  dat = phecodes.afr[, keepcols, with = FALSE]
  colnames(dat)[5] = "nvisits"
  dat[nvisits >= 2, `:=`(phenotype, 1)]
  dat[nvisits == 0, `:=`(phenotype, 0)]
  l1 = glm(data = dat, phenotype ~ Sex + age_cent + age_cent2, family = binomial(link = "logit"),
    control = list(maxit = 1000), na.action = "na.exclude")
  dat[, `:=`(residuals, resid(l1))]
  return(dat[, .(European, regional_haplo, residuals)])

})

names(bpheno.resids) = bphenotypes
bpheno.resids = dplyr::bind_rows(bpheno.resids, .id = "phenotype")

fwrite(bpheno.resids, "~/mitonuclear2/data/mtgwas/mitonuc.incomp/pmbb_bphenoresids_sex_age_age2_08112022.txt",
  sep = "\t")


keepcols = c("IID", "haplo_num", "Sex", "age_cent", "age_cent2", "African", "401")
ht = phecodes.afr[, keepcols, with = FALSE]
colnames(ht)[7] = "nvisits"
ht[nvisits >= 2, `:=`(phenotype, 1)]
ht[nvisits == 0, `:=`(phenotype, 0)]
ht = na.omit(ht)
l1 = glm(data = ht, phenotype ~ haplo_num + African + Sex + age_cent + age_cent2 +
  haplo_num * African, family = binomial(link = "logit"), control = list(maxit = 1000),
  na.action = "na.exclude")

pred.dat = ht[, .(haplo_num, African)]
pred.dat[, `:=`(Sex, "Female")]
pred.dat[, `:=`(age_cent, -0.03)]
pred.dat[, `:=`(age_cent2, 1)]
pred0 = predict(l1, pred.dat[haplo_num == 0, ], se.fit = TRUE, type = "response")
pred1 = predict(l1, pred.dat[haplo_num == 1, ], se.fit = TRUE, type = "response")

pred.dat[haplo_num == 0, `:=`(pred.y, pred0$fit)]
pred.dat[haplo_num == 1, `:=`(pred.y, pred1$fit)]

pred.dat[haplo_num == 0, `:=`(se.y, pred0$se.fit)]
pred.dat[haplo_num == 1, `:=`(se.y, pred1$se.fit)]

fwrite(pred.dat, "~/mitonuclear2/data/mtgwas/mitonuc.incomp/pmbb_htdat_08082022.txt",
  sep = "\t")


###### anemias

anemias = c("285", "285.2")
sapply(anemias, function(pheno.code) {
  keepcols = c("IID", "Sex", "age_cent", "age_cent2", pheno.code, "African", "regional_haplo")

  dat = phecodes.afr[, keepcols, with = FALSE]
  colnames(dat)[5] = "nvisits"
  dat[nvisits >= 2, `:=`(phenotype, 1)]
  dat[nvisits == 0, `:=`(phenotype, 0)]
  l1 = glm(data = dat, phenotype ~ Sex + age_cent + age_cent2, family = binomial(link = "logit"),
    control = list(maxit = 1000), na.action = "na.exclude")

  ht.dat = data.table(resid = resid(l1), African = dat$African, regional_haplo = dat$regional_haplo)

  fwrite(ht.dat, paste("~/mitonuclear2/data/mtgwas/mitonuc.incomp/pmbb_residuals_",
    pheno.code, "_sex_age_age2_08082022.txt", sep = ""), sep = "\t")

})

bphenos = c("070.2", "401.1", "401", "580.4","580.11","580")
bpheno.preds = lapply(bphenos, function(x) {
  keepcols = c("IID", "regional_haplo", "Sex", "age_cent", "age_cent2", "European",
    x)
  dat = phecodes.afr[, .SD, .SDcols = keepcols]
  colnames(dat)[7] = "nvisits"
  dat[nvisits >= 2, `:=`(phenotype, 1)]
  dat[nvisits == 0, `:=`(phenotype, 0)]
  dat = na.omit(dat)
  l1 = glm(data = dat, phenotype ~ regional_haplo + European + Sex + age_cent +
    age_cent2 + regional_haplo * European, family = binomial(link = "logit"),
    control = list(maxit = 1000), na.action = "na.exclude")

  pred.dat = dat[, .(regional_haplo, European)]
  pred.dat[, `:=`(Sex, "Female")]
  pred.dat[, `:=`(age_cent, -0.03)]
  pred.dat[, `:=`(age_cent2, 1)]
  pred0 = predict(l1, pred.dat[regional_haplo == "African", ], se.fit = TRUE, type = "response")
  pred1 = predict(l1, pred.dat[regional_haplo == "European", ], se.fit = TRUE,
    type = "response")

  pred.dat[regional_haplo == "African", `:=`(pred.y, pred0$fit)]
  pred.dat[regional_haplo == "European", `:=`(pred.y, pred1$fit)]

  pred.dat[regional_haplo == "African", `:=`(se.y, pred0$se.fit)]
  pred.dat[regional_haplo == "European", `:=`(se.y, pred1$se.fit)]

  return(pred.dat)

})

names(bpheno.preds) = bphenos
bpheno.preds = dplyr::bind_rows(bpheno.preds, .id = "Phenotype")

fwrite(bpheno.preds, "~/mitonuclear2/data/mtgwas/mitonuc.incomp/pmbb_bphenopreds_dat_08152022.txt",
  sep = "\t")
