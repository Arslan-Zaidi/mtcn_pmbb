


suppressWarnings(suppressMessages({
  library(data.table)
  library(pbapply)
}))

# define root directory F = rprojroot::is_rstudio_project$make_fix_file()

"%ni%" <- Negate("%in%")


# read ancestry information
ancestry = fread("~/mitonuclear2/data/genotype/pmbb.afr.weighted.04222022.glanc")


################# Quant phenotypes ###############



afr.covar = fread("~/mitonuclear2/data/phenotype/PMBB_plink4mat_AFR_07152022.covar")
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
afr.covar = merge(afr.covar, ancestry[, .(IID, African, European)], by = "IID")

# read median of the labs
labs = fread("~/mitonuclear2/data/phenotype/PMBB_2.2_quant_08122022.boxcox.pheno",
  header = TRUE)

  labs = merge(labs, afr.covar, by.x = "PMBB_ID", by.y = "IID")

# labs = labs[,replace(.SD, is.na(.SD)==TRUE, -9)]
# colnames(labs)[1] = "IID"
# labs = labs[IID %in% afr.covar$IID, ]
# labs = as.data.table(reshape2::melt(labs, id.vars = "IID"))
# labs = labs[is.na(value) == FALSE, ]
# labs2 = labs[, `:=`(int.value, qnorm((rank(value) - 0.5)/length(value))), by = "variable"]
# labs2 = labs2[variable %ni% c("Bands_percent", "Bands_THOperuL", "Basophils_percent",
#   "Basophils_THOperuL", "Eosinophil_percent", "Eosinophil_THOperuL", "Lymphocyte_percent",
#   "Lymphocyte_THOperuL", " Metamyeloctye", "Monocyte_percent", "Monocyte_THOperuL",
#   "Myelocyte", "PLT", "Promyelocyte", "WBC"), ]
#
# labs2 = dcast(labs2, IID ~ variable, value.var = "int.value")
# labs2 = merge(labs2, afr.covar, by = "IID")
#
# # process mtcn data - convert -9 to NA
# mtcn = fread("~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_02092022.txt")
# mtcn[mean.mt == -9, `:=`(mean.mt, NA)]
#
# labs2 = merge(labs2, mtcn, by = c("FID", "IID"))
# labs2[, `:=`(mean.mt, qnorm((rank(mean.mt) - 0.5)/length(mean.mt)))]


assumed.effects = seq(0.05, 1.5, 0.05)
# function to run GWAS
cal.pwr = function(phenotype, nboots) {

  dat = labs[Phenotype == phenotype,]
  colnames(dat)[8] = "phenotype"
  dat = na.omit(dat)
  l1 = lm(data = dat, tvalue ~ Sex + age_cent + age_cent2)

  s1 = summary(l1)
  ve = s1$sigma^2

  #dat[haplo_num == 0, `:=`(haplo_num, -1)]
  pwr = rep(NA, length(assumed.effects))
  ninds = nrow(dat)
  for (i in 1:length(assumed.effects)) {
    pvec = rep(NA, nboots)
    for (j in 1:nboots) {
      dat[regional_haplo == "European", sim.phenotype := assumed.effects[i] * African]
        dat[regional_haplo == "African", sim.phenotype := assumed.effects[i] * European]
        dat[, sim.phenotype := rnorm(ninds, 0, sqrt(ve))]
      l1 = lm(data = dat, sim.phenotype ~ European + regional_haplo + regional_haplo *
        European)
      s1 = summary(l1)$coefficients
      pvec[j] = s1[4, 4]
    }
    pwr[i] = length(which(pvec < 5e-04))/nboots

  }

  return(data.table(power = pwr))
}

lab.phenotypes = unique(labs$Phenotype)

pwrs = pblapply(lab.phenotypes, cal.pwr, nboots = 1e3)
names(pwrs) = lab.phenotypes
pwrs.df = dplyr::bind_rows(pwrs, .id = "phenotype")
pwrs.df$effect = assumed.effects

# # function to run GWAS
# cal.se = function(phenotype) {
#
#   keepcols = c("IID", "haplo_num", "Sex", "age_cent", "age_cent2", "African", phenotype)
#
#   dat = labs2[, keepcols, with = FALSE]
#   colnames(dat)[7] = "phenotype"
#   dat = na.omit(dat)
#   l1 = lm(data = dat, phenotype ~ haplo_num + African + Sex + age_cent + age_cent2 +
#     haplo_num * African)
#
#   s1 = summary(l1)
#   ve = s1$sigma^2
#
#   se = s1$coefficients[4, 2]
#   return(c(se, nrow(dat), ve))
# }


fwrite(pwrs.df, "~/mitonuclear2/data/mtgwas/mitonuc.incomp/pwrsim_incompatibility_labs_08302022.txt",
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# ses = pblapply(lab.phenotypes, cal.se)
# names(ses) = lab.phenotypes
# ses.df = bind_rows(ses, .id = "phenotype")
#
#
# fwrite(ses.df, "~/mitonuclear2/data/mtgwas/mitonuc.incomp/pwrses_labs.txt", sep = "\t",
#   col.names = TRUE, row.names = FALSE, quote = FALSE)
