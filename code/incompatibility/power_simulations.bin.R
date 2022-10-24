
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-p", "--phecode"), type = "character", help = "phecode"),
  make_option(c("-n","--nreplicates"), type = "integer", default = 1000, help = "no. of replicates")
)

args <- parse_args(OptionParser(option_list = option_list))

suppressWarnings(suppressMessages({
  library(data.table)
  library(pbapply)
}))

# define root directory F = rprojroot::is_rstudio_project$make_fix_file()

"%ni%" <- Negate("%in%")


# read ancestry information
ancestry = fread("~/mitonuclear2/data/genotype/pmbb.afr.weighted.04222022.glanc")


################# Quant phenotypes ###############



phecodes = fread("~/pmbb/Phenotype/2.2/PMBB-Release-2020-2.2_phenotype_PheCode-matrix.txt",
  header = TRUE)


bin.phenotypes = setdiff(colnames(phecodes), c("PMBB_ID"))

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

phecodes.afr = merge(afr.covar, phecodes, by.y = "PMBB_ID", by.x = "IID")

bin.phenotypes = setdiff(colnames(phecodes), c("PMBB_ID"))
phecodes.afr = merge(afr.covar, phecodes, by.y = "PMBB_ID", by.x = c("IID"))

#### read previous results

assumed.or = c(1.5, 2, 2.3, 3, 3.5, 4) # effect of African ancestry
neffects = length(assumed.or)

# function to run GWAS
cal.pwr.logistic = function(pheno.code, nboots) {

  keepcols = c("IID", "regional_haplo", "Sex", "age_cent", "age_cent2", "European", "African",pheno.code)

  dat = phecodes.afr[, keepcols, with = FALSE]
  colnames(dat)[8] = "case.control"
  dat[case.control == 0, `:=`(phenotype, 0)]
  dat[case.control > 2, `:=`(phenotype, 1)]
  l1 = glm(data = dat, phenotype ~ Sex + age_cent + age_cent2, family = binomial(link = "logit"),
    control = list(maxit = 1000))

  #dat[haplo_num == 0, `:=`(haplo_num, -1)]
  dat2 = na.omit(dat)

  # prevalence
  mean.sex = mean(dat$Sex == "Male")
  mean.age = mean(dat$age_cent)
  mean.age2 = mean(dat$age_cent2)

  b0 = sum(coefficients(l1) * c(1, mean.sex, mean.age, mean.age2))

  pwr = rep(NA, neffects)
  for (i in 1:neffects) {
    pvec = rep(NA, nboots)
    for (j in 1:nboots) {

      b1 = log(assumed.or[i])
      dat2[regional_haplo == "European", y0 := b0 + b1 * African]
      dat2[regional_haplo == "African", y0 := b0 - b1 * African]
      dat2[, probs := exp(y0)/(1 + exp(y0))]
      y = rbinom(nrow(dat2), 1, dat2$probs)
      dat2[, sim.phenotype := y]

      l2 = glm(data = dat2, sim.phenotype ~ European + regional_haplo + regional_haplo *
        European, family = "binomial")
      s2 = summary(l2)$coefficients
      pvec[j] = s2[4, 4]
    }
    pwr[i] = mean(pvec < 5e-05)

  }
  return(data.table(phecode = pheno.code, esize.or = assumed.or, power = pwr))
}

pwrs.df = cal.pwr.logistic(args$phecode, nboots = args$nreplicates)
# names(pwrs) = args$phecodec
# pwrs.df = dplyr::bind_rows(pwrs, .id = "phenotype")

output = paste("~/mitonuclear2/data/mtgwas/mitonuc.incomp/power/pwrsim_incompatibility_bin.", args$phecode, ".txt", sep = "")

fwrite(pwrs.df, output,
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
