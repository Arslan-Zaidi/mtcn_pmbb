library(data.table)
library(ggplot2)
library(pbapply)

"%ni%" <- Negate("%in%")



# read mtcn
mtcn = fread("~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_10212022.pheno")

# #read phecode data phecodes =
# fread('~/mitonuclear2/data/phenotype/PMBB_Diagnosis_Deidentified_05102022_Phecodes.txt.gz',
# header = TRUE)

phecodes = fread("~/pmbb/Phenotype/2.2/PMBB-Release-2020-2.2_phenotype_PheCode-matrix.txt",
  header = TRUE)
# phecodes.afr =
# fread('~/mitonuclear2/data/phenotype/PMBB_plink4mat_phecodes_AFR.pheno',
# header = TRUE) phecodes.afr = melt(phecodes.afr, id.vars = c('FID','IID'))
# phecodes.afr[value == -9, value := NA] phecodes.afr = dcast(phecodes.afr, FID
# + IID ~ variable)

bin.phenotypes = setdiff(colnames(phecodes), c("PMBB_ID"))

# read genetic PCs computed separately on each cohort
pcs.afr = fread("~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR.maf0.01.pca.eigenvec")
pcs.eur = fread("~/mitonuclear2/data/gcta/EUR/PMBB-Release-2020-2.0_genetic_genotype.EUR.maf0.01.pca.eigenvec")
colnames(pcs.afr) = colnames(pcs.eur) = c("FID", "IID", paste("PC", seq(1, 20), sep = ""))

# read sex and age information
demog = fread("~/pmbb/Phenotype/2.2/PMBB-Release-2020-2.2_phenotype_covariates.txt",
  header = TRUE)
# calculate age assuming collected in April 2020 demog[, Age := interval(start
# = ymd(BD), end = ymd('2020-04-01'))/duration(n=1, unit = 'years')]
# demog[GENDER_CODE == 'M', Sex := 'M'] demog[GENDER_CODE == 'F', Sex := 'F']
demog[is.na(Age_at_Enrollment) == FALSE, `:=`(Age, Age_at_Enrollment)]
demog[is.na(Age_at_Enrollment) == TRUE, `:=`(Age, Age_first_sample)]
demog = demog[Age > 0]

covar = demog[, .(PMBB_ID, Sex, Age)]

mtcn.afr = merge(mtcn, pcs.afr, by = c("FID", "IID"))
mtcn.eur = merge(mtcn, pcs.eur, by = c("FID", "IID"))

mtcn.afr = merge(mtcn.afr, covar, by.x = "IID", by.y = "PMBB_ID")
mtcn.eur = merge(mtcn.eur, covar, by.x = "IID", by.y = "PMBB_ID")

# merge all files
all.afr = merge(mtcn.afr, phecodes, by.x = "IID", by.y = "PMBB_ID")
all.eur = merge(mtcn.eur, phecodes, by.x = "IID", by.y = "PMBB_ID")


"%ni%" <- Negate("%in%")
pc.covars = paste("PC", seq(1, 20), sep = "")
pc.form = paste(pc.covars, collapse = "+")

phewas.glm = function(phecode, Pop) {
  if (Pop == "AFR") {
    dat = all.afr[, .SD, .SDcols = c("IID", phecode, "rlrmtcn", "Sex", "Age",
      pc.covars)]
  }
  if (Pop == "EUR") {
    dat = all.eur[, .SD, .SDcols = c("IID", phecode, "rlrmtcn", "Sex", "Age",
      pc.covars)]
  }
  colnames(dat)[2] = "nvisits"
  dat = na.omit(dat)

  dat[nvisits >= 2, `:=`(phenotype, 1)]
  dat[nvisits == 0, `:=`(phenotype, 0)]

  ncases = length(which(dat$phenotype == 1))
  ncontrols = length(which(dat$phenotype == 0))
  form1 = paste("phenotype ~ rlrmtcn + Sex + poly(Age, 2)+", pc.form)
  l1 = glm(data = dat, formula = as.formula(form1), family = binomial(link = "logit"))

  s1 = summary(l1)$coefficients

  lm1 = data.table(Beta = s1[, 1], SE = s1[, 2], Z = s1[, 3], P = s1[, 4], ncases = ncases,
    ncontrols = ncontrols)
  lm1$Variable = rownames(s1)

  lm1 = lm1[Variable %ni% paste("PC", seq(1, 20), sep = ""), ]

  lm1$Pop = Pop
  return(lm1)
}


glm.afr.list = pblapply(bin.phenotypes, function(x) {
  phewas.glm(x, Pop = "AFR")
})
names(glm.afr.list) = bin.phenotypes
glm.afr.df = dplyr::bind_rows(glm.afr.list, .id = "Phecode")
fwrite(glm.afr.df, "~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.phewas_phecodes.afr.10212022.glm.logistic")

glm.eur.list = pblapply(bin.phenotypes, function(x) {
  phewas.glm(x, Pop = "EUR")
})
names(glm.eur.list) = bin.phenotypes
glm.eur.df = dplyr::bind_rows(glm.eur.list, .id = "Phecode")
fwrite(glm.eur.df, "~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.phewas_phecodes.eur.10212022.glm.logistic")



############## associations with labs (INT) ####################### labs.afr =
############## fread('~/mitonuclear2/data/phenotype/PMBB_plink4mat_labs_INT_05102022.AFR.pheno')
############## labs.eur =
############## fread('~/mitonuclear2/data/phenotype/PMBB_plink4mat_labs_INT_05102022.EUR.pheno')
############## lab.cols = c('A1C', 'ALB', 'Alkphos', 'ALT', 'AST',
############## 'Bili_Direct', 'Bili_Indirect', 'Bili_Total', 'CK', 'CRP',
############## 'Chol_HDL', 'Chol_LDL', 'Chol_Total', 'Creatinine', 'ESR',
############## 'GGT', 'Glucose', 'Glucose_fasting', 'INR', 'PTT',
############## 'triglycerides', 'urobili') labs.afr = labs.afr[, .SD, .SDcols =
############## c('FID', 'IID', lab.cols)] labs.afr = melt(labs.afr, id.vars =
############## c('FID', 'IID')) labs.afr[value == -9, `:=`(value, NA)] labs.afr
############## = dcast(labs.afr, FID + IID ~ variable) labs.afr =
############## merge(labs.afr, mtcn.afr, by = c('FID', 'IID')) labs.eur =
############## labs.eur[, .SD, .SDcols = c('FID', 'IID', lab.cols)] labs.eur =
############## melt(labs.eur, id.vars = c('FID', 'IID')) labs.eur[value == -9,
############## `:=`(value, NA)] labs.eur = dcast(labs.eur, FID + IID ~
############## variable) labs.eur = merge(labs.eur, mtcn.eur, by = c('FID',
############## 'IID')) # wbc.afr = #
############## fread('~/mitonuclear2/data/phenotype/PMBB_plink4mat_WBC_INT_05102022.AFR.pheno')
############## # wbc.eur = #
############## fread('~/mitonuclear2/data/phenotype/PMBB_plink4mat_WBC_INT_05102022.EUR.pheno')
############## lab.cols = setdiff(lab.cols, c('FID', 'IID')) # labs.afr =
############## merge(labs.afr, covar, by = c('FID','IID')) labs.eur = #
############## merge(labs.eur, covar, by = c('FID','IID'))

# lm.afr.list = pblapply(quant.cols, phewas.labs, Pop = 'AFR') lm.eur.list =
# pblapply(quant.cols, phewas.labs, Pop = 'EUR') names(lm.afr.list) =
# names(lm.eur.list) = quant.cols lm.afr = dplyr::bind_rows(lm.afr.list, .id =
# 'Phenotype') lm.eur = dplyr::bind_rows(lm.eur.list, .id = 'Phenotype')
# fwrite(lm.afr,
# '~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.phewas_labs.int.afr.08042022.linear')
# fwrite(lm.eur,
# '~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.phewas_labs.int.eur.08042022.linear')


####### association with labs box-cox transformed #####

quants = fread("~/mitonuclear2/data/phenotype/PMBB_2.2_quant_08122022.boxcox.pheno",
  header = TRUE)
quants = quants[Phenotype %ni% c("Basophil_median", "Eosinophil_median", "Lymophocyte_median",
  "Monocyte_median", "Neutrophil_median", "PLATELET_median", "WBC_median")]

labs.eur = merge(mtcn.eur, quants, by.x = "IID", by.y = "PMBB_ID")
labs.afr = merge(mtcn.afr, quants, by.x = "IID", by.y = "PMBB_ID")

labs.eur[, `:=`(scaled_value, scale(tvalue)), by = "Phenotype"]
labs.afr[, `:=`(scaled_value, scale(tvalue)), by = "Phenotype"]
quant.cols = unique(quants$Phenotype)


phewas.labs = function(lab, Pop) {
  if (Pop == "AFR") {
    dat = labs.afr[Phenotype == lab]
  }
  if (Pop == "EUR") {
    dat = labs.eur[Phenotype == lab]
  }
  dat = na.omit(dat)
  nsamples = nrow(dat)

  form1 = paste("scaled_value ~ rlrmtcn+ Sex + poly(Age, 2)+", pc.form)
  l1 = lm(data = dat, formula = as.formula(form1))

  s1 = summary(l1)$coefficients

  lm1 = data.table(Beta = s1[, 1], SE = s1[, 2], Z = s1[, 3], P = s1[, 4], nsamples = nsamples)
  lm1$Variable = rownames(s1)

  lm1 = lm1[Variable %ni% paste("PC", seq(1, 20), sep = ""), ]

  lm1$Pop = Pop
  return(lm1)

}

lm.afr.list = pblapply(quant.cols, phewas.labs, Pop = "AFR")
lm.eur.list = pblapply(quant.cols, phewas.labs, Pop = "EUR")

names(lm.afr.list) = names(lm.eur.list) = quant.cols
lm.afr = dplyr::bind_rows(lm.afr.list, .id = "Phenotype")
lm.eur = dplyr::bind_rows(lm.eur.list, .id = "Phenotype")

fwrite(lm.afr, "~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.phewas_quant.bxcx.afr.10212022.linear")
fwrite(lm.eur, "~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.phewas_quant.bxcx.eur.10212022.linear")

#### INT labs.eur[, `:=`(int.value, qnorm((rank(value) - 0.5)/length(value))),
#### by = 'variable'] labs.afr[, `:=`(int.value, qnorm((rank(value) -
#### 0.5)/length(value))), by = 'variable'] phewas.labs.int = function(lab,
#### Pop) { if (Pop == 'AFR') { dat = labs.afr[variable == lab] } if (Pop ==
#### 'EUR') { dat = labs.eur[variable == lab] } dat = na.omit(dat) nsamples =
#### nrow(dat) form1 = paste('int.value ~ rlrmtcn+ Sex + poly(Age, 2, raw =
#### TRUE)+', pc.form) l1 = lm(data = dat, formula = as.formula(form1)) s1 =
#### summary(l1)$coefficients lm1 = data.table(Beta = s1[, 1], SE = s1[, 2], Z
#### = s1[, 3], P = s1[, 4], nsamples = nsamples) lm1$Variable = rownames(s1)
#### lm1 = lm1[Variable %ni% paste('PC', seq(1, 20), sep = ''), ] lm1$Pop = Pop
#### return(lm1) } lm.afr.list = pblapply(quant.cols, phewas.labs.int, Pop =
#### 'AFR') lm.eur.list = pblapply(quant.cols, phewas.labs.int, Pop = 'EUR')
#### names(lm.afr.list) = names(lm.eur.list) = quant.cols lm.afr =
#### dplyr::bind_rows(lm.afr.list, .id = 'Phenotype') lm.eur =
#### dplyr::bind_rows(lm.eur.list, .id = 'Phenotype') fwrite(lm.afr,
#### '~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.phewas_quant.int.afr.08112022.linear')
#### fwrite(lm.eur,
#### '~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.phewas_quant.int.eur.08112022.linear')


########################## using untransformed values
########################## ########################## mtcn =
########################## fread('~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_mean.mtcn.residuals_08042022.pheno')
########################## colnames(mtcn)[3] = 'rlrmtcn' # just to make the
########################## script work mtcn.afr = merge(mtcn, pcs.afr, by =
########################## c('FID', 'IID')) mtcn.eur = merge(mtcn, pcs.eur, by
########################## = c('FID', 'IID')) mtcn.afr = merge(mtcn.afr, covar,
########################## by.x = 'IID', by.y = 'PMBB_ID') mtcn.eur =
########################## merge(mtcn.eur, covar, by.x = 'IID', by.y =
########################## 'PMBB_ID') # merge all files all.afr =
########################## merge(mtcn.afr, phecodes, by.x = 'IID', by.y =
########################## 'PMBB_ID') all.eur = merge(mtcn.eur, phecodes, by.x
########################## = 'IID', by.y = 'PMBB_ID') glm.afr.list =
########################## pblapply(bin.phenotypes, function(x) { phewas.glm(x,
########################## Pop = 'AFR') }) names(glm.afr.list) = phenotypes
########################## glm.afr.df = dplyr::bind_rows(glm.afr.list, .id =
########################## 'Phecode') fwrite(glm.afr.df,
########################## '~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.untransformed.phewas_phecodes.afr.08042022.glm.logistic')
########################## glm.eur.list = pblapply(bin.phenotypes, function(x)
########################## { phewas.glm(x, Pop = 'EUR') }) names(glm.eur.list)
########################## = phenotypes glm.eur.df =
########################## dplyr::bind_rows(glm.eur.list, .id = 'Phecode')
########################## fwrite(glm.eur.df,
########################## '~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_mtcn.untransformed.phewas_phecodes.eur.08042022.glm.logistic')
########################## labs =
########################## fread('~/mitonuclear2/data/phenotype/PMBB_labs_2.2_07282022.txt')
########################## labs = labs[LAB %in% lab.cols, ] labs[,
########################## `:=`(scaled_value, scale(MEDIAN_VALUE))] labs =
########################## labs[scaled_value < 10, ] labs = labs[,
########################## `:=`(scaled_value, scale(scaled_value))] labs =
########################## dcast(labs, PMBB_ID ~ LAB, value.var =
########################## 'scaled_value') labs.eur = merge(mtcn.eur, labs,
########################## by.x = 'IID', by.y = 'PMBB_ID') labs.afr =
########################## merge(mtcn.afr, labs, by.x = 'IID', by.y =
########################## 'PMBB_ID') lm.afr.list = pblapply(lab.cols,
########################## phewas.labs, Pop = 'AFR') lm.eur.list =
########################## pblapply(lab.cols, phewas.labs, Pop = 'EUR')
########################## names(lm.afr.list) = names(lm.eur.list) = lab.cols
########################## lm.afr = dplyr::bind_rows(lm.afr.list, .id = 'Lab')
########################## lm.eur = dplyr::bind_rows(lm.eur.list, .id = 'Lab')
########################## fwrite(lm.afr,
########################## '~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_untransformed.mtcn.phewas_labs.raw.afr.08042022.linear')
########################## fwrite(lm.eur,
########################## '~/mitonuclear2/data/mtgwas/mtcn.phewas/pmbb_untransformed.mtcn.phewas_labs.raw.eur.08042022.linear')


### follow up with cardiovascular traits #####

# using Hagg et al model to test for association between mtDNA copy number and
# cardiovascular traits
cvd.traits = c("394", "394.2", "395", "395.2", "411", "425", "427", "427.1", "427.2",
  "427.21", "427.22", "427.6", "427.9")
mall.eur = melt(all.eur[, .SD, .SDcols = c("IID", "FID", "rlrmtcn", "Age", "Sex",
  pc.covars, cvd.traits)], id.vars = c("IID", "FID", "rlrmtcn", "Age", "Sex", pc.covars))
colnames(mall.eur)[27] = "nvisits"
mall.eur[nvisits >= 2, `:=`(phenotype, 1)]
mall.eur[nvisits == 0, `:=`(phenotype, 0)]

wbc = fread("~/mitonuclear2/data/phenotype/WBC_summary_data_with_platelets.csv")
wbc[, `:=`(Neutrophil_p, Neutrophil_median/WBC_median)]
wbc[, `:=`(Lymphocyte_p, Lymphocyte_median/WBC_median)]

mall.eur = merge(mall.eur, wbc[, .(PMBB_ID, Neutrophil_p, Lymphocyte_p, WBC_median,
  PLATELET_median)], by.x = "IID", by.y = "PMBB_ID")
mtcn = fread("~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_07122022.txt")
mall.eur = merge(mall.eur, mtcn, by = c("FID", "IID"))
l1 = lapply(cvd.traits, function(x) {
  dat = mall.eur[variable == x]
  form1 = paste("phenotype ~ log.mtcn + Sex + Age+ Neutrophil_p + Lymphocyte_p + WBC_median")
  l1 = glm(data = dat, form1, family = binomial(link = "logit"))
  s1 = summary(l1)$coefficients
  lm1 = data.table(Beta = s1[, 1], SE = s1[, 2], Z = s1[, 3], P = s1[, 4])
  lm1$Variable = rownames(s1)
  lm1 = lm1[Variable %ni% paste("PC", seq(1, 20), sep = ""), ]
  lm1$Model = "Hagg"

  form2 = paste("phenotype ~ log.mtcn + Sex + Age+ Neutrophil_p + Lymphocyte_p + WBC_median+ PLATELET_median")
  l2 = glm(data = dat, form2, family = binomial(link = "logit"))
  s2 = summary(l2)$coefficients
  lm2 = data.table(Beta = s2[, 1], SE = s2[, 2], Z = s2[, 3], P = s2[, 4])
  lm2$Variable = rownames(s2)
  lm2 = lm2[Variable %ni% paste("PC", seq(1, 20), sep = ""), ]
  lm2$Model = "Hagg + PLTs"

  lm3 = rbind(lm1, lm2)

  return(lm3)
})


names(l1) = cvd.traits
l1 = dplyr::bind_rows(l1, .id = "Phecode")
l1[Variable == "log.mtcn"]


##### follow up with liver traits
liver.damage = c("571.8", "571.81", "530.2", "573", "571", "573.2", "571.51")
alcoholism = c("317", "317.1", "317.11")

# check if alcoholism is basically the same
phecodes[, lapply(.SD, mean), .SDcols = alcoholism]
# 317 317.1 317.11 1: 0.09301847 0.09272099 0.09219469

# just use one in the correction
liver.associations = lapply(liver.damage, function(x) {
  dat = all.eur[, .SD, .SDcols = c("IID", x, "rlrmtcn", "317", "Sex", "Age", pc.covars)]
  colnames(dat)[c(2, 4)] = c("nvisits", "alcoholism")
  dat = na.omit(dat)

  dat[nvisits >= 2, `:=`(phenotype, 1)]
  dat[nvisits == 0, `:=`(phenotype, 0)]
  dat[alcoholism >= 2, `:=`(covariate, 1)]
  dat[alcoholism == 0, `:=`(covariate, 0)]

  form1 = paste("phenotype ~ rlrmtcn + Sex + poly(Age, 2)+", pc.form)
  l1 = glm(data = dat, formula = as.formula(form1), family = binomial(link = "logit"))
  s1 = summary(l1)$coefficients
  lm1 = data.table(Beta = s1[, 1], SE = s1[, 2], Z = s1[, 3], P = s1[, 4])
  lm1$Variable = rownames(s1)
  lm1 = lm1[Variable %ni% paste("PC", seq(1, 20), sep = ""), ]
  lm1$Model = "no alocholism"

  form2 = paste("phenotype ~ rlrmtcn + as.factor(covariate) + Sex + poly(Age, 2)+",
    pc.form)
  l2 = glm(data = dat, formula = as.formula(form2), family = binomial(link = "logit"))
  s2 = summary(l2)$coefficients
  lm2 = data.table(Beta = s2[, 1], SE = s2[, 2], Z = s2[, 3], P = s2[, 4])
  lm2$Variable = rownames(s2)
  lm2 = lm2[Variable %ni% paste("PC", seq(1, 20), sep = ""), ]
  lm2$Model = "alcoholism"

  lm3 = rbind(lm1, lm2)
  return(lm3)
})

names(liver.associations) = liver.damage
liver.associations = dplyr::bind_rows(liver.associations, .id = "Phecode")
library(PheWAS)
liver.associations = addPhecodeInfo(liver.associations)
liver.associations = liver.associations[Variable == "rlrmtcn"]

fwrite(liver.associations, "~/mitonuclear2/output/tables/supptable_mtcn_phewas_liver_10212022.txt",
  sep = "\t")
