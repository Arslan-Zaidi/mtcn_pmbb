
library(data.table)


### covariate file
demog = fread("~/pmbb/Phenotype/2.2/PMBB-Release-2020-2.2_phenotype_covariates.txt",
  header = TRUE)
# calculate age assuming collected in April 2020 demog[, Age := interval(start
# = ymd(BD), end = ymd('2020-04-01'))/duration(n=1, unit = 'years')]
# demog[GENDER_CODE == 'M', Sex := 'M'] demog[GENDER_CODE == 'F', Sex := 'F']
demog[is.na(Age_at_Enrollment) == FALSE, `:=`(Age, Age_at_Enrollment)]
demog[is.na(Age_at_Enrollment) == TRUE, `:=`(Age, Age_first_sample)]
demog = demog[Age > 0]

covar = demog[, .(PMBB_ID, Sex, Age)]
# covar = demog[,Age2 := Age^2]


# genetic PCs 1 - 20h
afr.pca = fread("~/mitonuclear2/data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR.maf0.01.pca.eigenvec")
afr.pca = afr.pca[, c(2:22)]
colnames(afr.pca) = c("IID", paste("PC", seq(1, 20), sep = ""))
afr.covar = merge(covar, afr.pca, by.x = "PMBB_ID", by.y = "IID")
#get mean and sd of age in thge AFR cohort (for later)
with(afr.covar, mean(Age, na.rm = TRUE)) #51.97
with(afr.covar, sd(Age, na.rm = TRUE)) #15.91
afr.covar[, `:=`(age_cent, (Age - mean(Age, na.rm = TRUE))/sd(Age, na.rm = TRUE))]  #standardize age
afr.covar[, `:=`(age_cent2, age_cent^2)]
afr.covar$Age = NULL
afr.covar = cbind(afr.covar$PMBB_ID, afr.covar)

eur.pca = fread("~/mitonuclear2/data/gcta/EUR/PMBB-Release-2020-2.0_genetic_genotype.EUR.maf0.01.pca.eigenvec")
eur.pca = eur.pca[, c(2:22)]
colnames(eur.pca) = c("IID", paste("PC", seq(1, 20), sep = ""))
eur.covar = merge(covar, eur.pca, by.x = "PMBB_ID", by.y = "IID")
with(eur.covar, mean(Age, na.rm = TRUE)) #57.42
with(eur.covar, sd(Age, na.rm = TRUE)) #15.96
eur.covar[, `:=`(age_cent, (Age - mean(Age, na.rm = TRUE))/sd(Age, na.rm = TRUE))]  #standardize age
eur.covar[, `:=`(age_cent2, age_cent^2)]
eur.covar$Age = NULL
eur.covar = cbind(eur.covar$PMBB_ID, eur.covar)

colnames(afr.covar)[c(1, 2)] = colnames(eur.covar)[c(1, 2)] = c("FID", "IID")

haplo = fread("~/mitonuclear2/data/genotype/PMBB_Release-2020-2.0-MT.haplo.gz")
colnames(haplo)[1] = c("IID")
haplo = haplo[, c(1:2)]
haplo$FID = haplo$IID
haplo = haplo[, .(FID, IID, Haplogroup)]
haplo$Haplo1 = substr(haplo$Haplogroup, 1, 1)
haplo$Haplo2 = substr(haplo$Haplogroup, 1, 2)

haplo.eur = haplo[IID %in% eur.covar$IID]
haplo.afr = haplo[IID %in% afr.covar$IID]

afr.covar = merge(afr.covar, haplo.afr, by = c("FID", "IID"))
eur.covar = merge(eur.covar, haplo.eur, by = c("FID", "IID"))

fwrite(eur.covar, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_EUR_07152022.covar",
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

fwrite(afr.covar, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_AFR_07152022.covar",
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


### phenotype file
phecodes = fread("~/mitonuclear2/data/phenotype/PMBB_Diagnosis_Deidentified_05102022_Phecodes.txt.gz",
  header = TRUE)

# relabel cases as 2 and controls as 1
phecodes = phecodes[, replace(.SD, .SD == 1, 2)]
phecodes = phecodes[, replace(.SD, .SD == 0, 1)]
phecodes = phecodes[, replace(.SD, is.na(.SD) == TRUE, -9)]
phecodes = cbind(phecodes$PMBB_ID, phecodes)
colnames(phecodes)[c(1, 2)] = c("FID", "IID")
#
# phecodes.afr = phecodes[IID %in% (afr.covar$IID)]
# phecodes.eur = phecodes[IID %in% c(eur.covar$IID)]

# phecodes.afr.m = as.data.table(reshape2::melt(phecodes.afr, id.vars =
# c('FID','IID'))) phecodes.afr.summary = phecodes.afr.m[,.(ncases =
# length(which(value == 2)), ntotal = length(which(value != -9))), by =
# 'variable'] fwrite(phecodes.afr.summary,
# '~/mitonuclear2/data/phenotype/PMBB_phecodes.AFR.05102022.summary', col.names
# = TRUE, row.names = FALSE, quote = FALSE, sep = '\t') phecodes.afr.summary =
# phecodes.afr.summary[ncases > 20,] phecodes.afr =
# phecodes.afr[,c('FID','IID',as.character(phecodes.afr.summary$variable)),
# with = FALSE] phecodes.eur.m = as.data.table(reshape2::melt(phecodes.eur,
# id.vars = c('FID','IID'))) phecodes.eur.summary = phecodes.eur.m[,.(ncases =
# length(which(value == 2)), ntotal = length(value != -9)), by = 'variable']
# fwrite(phecodes.eur.summary,
# '~/mitonuclear2/data/phenotype/PMBB_phecodes.EUR.summary', col.names = TRUE,
# row.names = FALSE, quote = FALSE, sep = '\t') phecodes.eur.summary =
# phecodes.eur.summary[ncases>20,] phecodes.eur = phecodes.eur[,
# c('FID','IID',as.character(phecodes.eur.summary$variable)), with = FALSE]
# fwrite(phecodes.afr,
# '~/mitonuclear2/data/phenotype/PMBB_plink4mat_phecodes_AFR.pheno', col.names
# = TRUE, row.names = FALSE, quote = FALSE, sep = '\t') fwrite(phecodes.eur,
# '~/mitonuclear2/data/phenotype/PMBB_plink4mat_phecodes_EUR.pheno', col.names
# = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


#### labs
labfiles = list.files(path = "~/pmbb/Phenotype/2.2", pattern = "PMBB-Release-2020-2.2_phenotype_labs-")
labfiles = labfiles[-grep("procedures", labfiles)]
labfiles = labfiles[-grep("PER", labfiles)]

labdat.list = list()
for (i in 1:length(labfiles)) {
  lab.suffix = unlist(strsplit(labfiles[i], split = "-"))
  lab.suffix = lab.suffix[-c(1:4)]
  if (length(lab.suffix) == 1) {
    lab = unlist(strsplit(lab.suffix, split = "\\."))[1]
  } else {
    lab.suffix2 = unlist(strsplit(lab.suffix[length(lab.suffix)], split = "\\."))[1]
    lab = paste(c(lab.suffix[-length(lab.suffix)], lab.suffix2), collapse = "_")
  }
  labdat = fread(paste("~/pmbb/Phenotype/2.2/", labfiles[i], sep = ""), header = TRUE)
  labdat = labdat[, .(PMBB_ID, RESULT_VALUE_NUM)]
  labdat[, `:=`(LAB, lab)]
  labdat.list[[i]] = labdat[, .(MEDIAN_VALUE = median(RESULT_VALUE_NUM, na.rm = TRUE)),
    by = c("PMBB_ID", "LAB")]
}

labs = dplyr::bind_rows(labdat.list)
fwrite(labs, "~/mitonuclear2/data/phenotype/PMBB_labs_2.2_07282022.txt", sep = "\t")


# inverse-rank normal transformation (separately in AFR and EUR)
labs.afr = labs[PMBB_ID %in% afr.covar$IID, ]
labs.eur = labs[PMBB_ID %in% eur.covar$IID, ]

labs.afr[, `:=`(int.value, qnorm((rank(MEDIAN_VALUE) - 0.5)/length(MEDIAN_VALUE))),
  by = "LAB"]
labs.eur[, `:=`(int.value, qnorm((rank(MEDIAN_VALUE) - 0.5)/length(MEDIAN_VALUE))),
  by = "LAB"]

labs.afr = dcast(labs.afr, PMBB_ID ~ LAB, value.var = "int.value")
labs.eur = dcast(labs.eur, PMBB_ID ~ LAB, value.var = "int.value")

labs.afr = labs.afr[, replace(.SD, is.na(.SD) == TRUE, -9)]
labs.eur = labs.eur[, replace(.SD, is.na(.SD) == TRUE, -9)]

labs.afr = cbind(labs.afr$PMBB_ID, labs.afr)
labs.eur = cbind(labs.eur$PMBB_ID, labs.eur)

colnames(labs.afr)[c(1, 2)] = colnames(labs.eur)[c(1, 2)] = c("FID", "IID")
labs2 = rbind(labs.afr, labs.eur)

fwrite(labs2, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_labs_INT_05102022.pheno",
  col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

fwrite(labs.afr, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_labs_INT_05102022.AFR.pheno",
  col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

fwrite(labs.eur, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_labs_INT_05102022.EUR.pheno",
  col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
##### inverse rank normal transofmration complete


# labs =
# fread('~/pmbb/Phenotype/PMBB-Release-2020-2.0_phenotype_labs-median-matrix.txt',
# header = TRUE) #remove all blood measurements from this - to be loaded from a
# different file labs = labs[, -grep('percent', colnames(labs)), with = FALSE]
# labs = labs[, -grep('THOperuL', colnames(labs)), with = FALSE] #load blood
# counts wbc = fread('~/mitonuclear2/data/phenotype/WBC_summary_data.csv')
# #keep median measurements wbc = wbc[, -grep('average', colnames(wbc)), with =
# FALSE] wbc = wbc[, -grep('date', colnames(wbc)), with = FALSE] labs =
# merge(labs, wbc, by = 'PMBB_ID')

labs = dcast(labs, PMBB_ID ~ LAB, value.var = "MEDIAN_VALUE")
fwrite(labs, "~/mitonuclear2/data/phenotype/PMBB_Labs.NA.05102022.pheno", sep = "\t")

labs = cbind(labs$PMBB_ID, labs)
colnames(labs)[c(1, 2)] = c("FID", "IID")

labs = labs[, replace(.SD, is.na(.SD) == TRUE, -9)]

fwrite(labs, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_labs_05102022.pheno",
  sep = "\t")

#### vitals
bp = fread("~/pmbb/Phenotype/2.2/PMBB-Release-2020-2.2_phenotype_vitals-BP.txt")
height = fread("~/pmbb/Phenotype/2.2/PMBB-Release-2020-2.2_phenotype_vitals-height.txt")
bmi = fread("~/pmbb/Phenotype/2.2/PMBB-Release-2020-2.2_phenotype_vitals-BMI.txt")

bp = bp[, .(PMBB_ID, SYSTOLIC, DIASTOLIC)]
bp = melt(bp, id.vars = "PMBB_ID")
bp = bp[, .(value = median(value)), by = c("PMBB_ID", "variable")]

height = height[, .(PMBB_ID, HEIGHT_IN)]
height = height[, .(value = median(HEIGHT_IN)), by = "PMBB_ID"]
height[, variable := "Height"]
height = height[, .(PMBB_ID, variable, value)]

bmi = bmi[, .(PMBB_ID, BMI)]
bmi = bmi[, .(value = median(BMI)), by = c("PMBB_ID")]
bmi[, variable := "BMI"]
bmi = bmi[, .(PMBB_ID, variable, value)]

vitals = rbind(bp, height, bmi)

fwrite(vitals, "~/mitonuclear2/data/phenotype/PMBB_vitals_2.2_08112022.txt", sep = "\t")


####### make covariate file with blood measurements (for mtdna copy number
####### GWAS) mlabs = melt(labs, id.vars = c('FID','IID')) mlabs = mlabs[,
####### .(nfull = length(which(value!=-9))), by = 'variable']

# blood.cols =
# c('Basophil_CON','Eosinophil_CON','Lymphocyte_CON','Monocyte_CON','Platelets','WBC')

wbc = fread("~/mitonuclear2/data/phenotype/WBC_summary_data_with_platelets.csv")

blood.cols = c("Basophil_median", "Eosinophil_median", "Lymphocyte_median", "Monocyte_median",
  "Neutrophil_median", "PLATELET_median", "WBC_median")

wbc = wbc[, c("PMBB_ID", blood.cols), with = FALSE]

wbc.afr = wbc[PMBB_ID %in% afr.covar$IID, ]
wbc.eur = wbc[PMBB_ID %in% eur.covar$IID, ]

wbc.afr = melt(wbc.afr, id.vars = c("PMBB_ID"), variable.name = "LAB", value.name = "VALUE")
wbc.eur = melt(wbc.eur, id.vars = c("PMBB_ID"), variable.name = "LAB", value.name = "VALUE")

wbc.afr2 = wbc.afr[, `:=`(int.value, qnorm((rank(VALUE) - 0.5)/length(VALUE))), by = "LAB"]
wbc.eur2 = wbc.eur[, `:=`(int.value, qnorm((rank(VALUE) - 0.5)/length(VALUE))), by = "LAB"]

wbc.afr2 = dcast(wbc.afr2, PMBB_ID ~ LAB, value.var = "int.value")
wbc.eur2 = dcast(wbc.eur2, PMBB_ID ~ LAB, value.var = "int.value")

wbc.afr2 = wbc.afr2[, replace(.SD, is.na(.SD) == TRUE, -9)]
wbc.eur2 = wbc.eur2[, replace(.SD, is.na(.SD) == TRUE, -9)]

wbc.afr2 = cbind(wbc.afr2$PMBB_ID, wbc.afr2)
wbc.eur2 = cbind(wbc.eur2$PMBB_ID, wbc.eur2)

colnames(wbc.afr2)[c(1, 2)] = colnames(wbc.eur2)[c(1, 2)] = c("FID", "IID")

# wbc = wbc[, (blood.cols):=lapply(.SD, function(x) { ifelse(is.na(x), -9, (x -
# mean(x,na.rm=TRUE))/sd(x, na.rm=TRUE)) }), .SDcols = blood.cols]

fwrite(wbc.afr2, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_WBC_INT_05102022.AFR.pheno",
  sep = "\t")
fwrite(wbc.eur2, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_WBC_INT_05102022.EUR.pheno",
  sep = "\t")

#get means and sd for blood counts
wbc.afr.scales = wbc.afr[, .(mean = mean(VALUE, na.rm =TRUE), sd = sd(VALUE, na.rm = TRUE)), by = "LAB"]
# LAB         mean          sd
# 1:   Basophil_median   0.02514477  0.03268704
# 2: Eosinophil_median   0.14356704  0.12441314
# 3: Lymphocyte_median   1.82310761  0.81524777
# 4:   Monocyte_median   0.54976142  0.19497530
# 5: Neutrophil_median   4.28594510  1.92064856
# 6:   PLATELET_median 241.43532819 75.71349929
# 7:        WBC_median   7.25193756  2.60436651

wbc.eur.scales = wbc.eur[, .(mean = mean(VALUE, na.rm =TRUE), sd = sd(VALUE, na.rm = TRUE)), by = "LAB"]
# LAB         mean          sd
# 1:   Basophil_median   0.03116015  0.05305579
# 2: Eosinophil_median   0.15591989  0.15817070
# 3: Lymphocyte_median   1.65248496  2.34769060
# 4:   Monocyte_median   0.57372136  0.23561531
# 5: Neutrophil_median   4.71541038  2.00517676
# 6:   PLATELET_median 224.75225234 73.14115702
# 7:        WBC_median   7.58190842  3.89910867

fwrite(wbc.afr.scales, "~/mitonuclear2/data/phenotype/PMBB_WBC_07212022.AFR.scales", sep = "\t")
fwrite(wbc.eur.scales, "~/mitonuclear2/data/phenotype/PMBB_WBC_07212022.EUR.scales", sep = "\t")


wbc.afr[, `:=`(VALUE, (VALUE - mean(VALUE, na.rm = TRUE))/sd(VALUE, na.rm = TRUE)), by = "LAB"]
wbc.eur[, `:=`(VALUE, (VALUE - mean(VALUE, na.rm = TRUE))/sd(VALUE, na.rm = TRUE)), by = "LAB"]

wbc.afr = dcast(wbc.afr, PMBB_ID ~ LAB, value.var = "VALUE")
wbc.eur = dcast(wbc.eur, PMBB_ID ~ LAB, value.var = "VALUE")

wbc.afr = wbc.afr[, replace(.SD, is.na(.SD) == TRUE, -9)]
wbc.eur = wbc.eur[, replace(.SD, is.na(.SD) == TRUE, -9)]

wbc.afr = cbind(wbc.afr$PMBB_ID, wbc.afr)
wbc.eur = cbind(wbc.eur$PMBB_ID, wbc.eur)

colnames(wbc.afr)[c(1, 2)] = colnames(wbc.eur)[c(1, 2)] = c("FID", "IID")

fwrite(wbc.afr, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_WBC_07152022.AFR.pheno",
  sep = "\t")
fwrite(wbc.eur, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_WBC_07152022.EUR.pheno",
  sep = "\t")


afr.covar2 = merge(afr.covar, wbc.afr, by = c("FID", "IID"))
eur.covar2 = merge(eur.covar, wbc.eur, by = c("FID", "IID"))

fwrite(eur.covar2, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_EUR.blood_07152022.covar",
  sep = "\t")
fwrite(afr.covar2, "~/mitonuclear2/data/phenotype/PMBB_plink4mat_AFR.blood_07152022.covar",
  sep = "\t")

# #print sample size for each lab measurement mlabs.afr = melt(labs.afr,
# id.vars = c('FID', 'IID')) mlabs.eur = melt(labs.eur, id.vars = c('FID',
# 'IID')) mlabs.afr.summary = mlabs.afr[, .(n = length(which(value != -9))), by
# = 'variable'] mlabs.eur.summary = mlabs.eur[, .(n = length(which(value !=
# -9))), by = 'variable'] fwrite(mlabs.afr.summary,
# '~/mitonuclear2/data/phenotype/PMBB_labs_03022022.n.AFR.txt', sep = '\t')
# fwrite(mlabs.eur.summary,
# '~/mitonuclear2/data/phenotype/PMBB_labs_03022022.n.EUR.txt', sep = '\t')

haplo = fread("~/mitonuclear2/data/genotype/PMBB_Release-2020-2.0-MT.haplo")
colnames(haplo)[1] = c("IID")
haplo = haplo[, c(1:2)]
haplo$FID = haplo$IID
haplo = haplo[, .(FID, IID, Haplogroup)]
haplo$Haplo1 = substr(haplo$Haplogroup, 1, 1)
haplo$Haplo2 = substr(haplo$Haplogroup, 1, 2)

haplo.eur = haplo[IID %in% eur.covar$IID]
haplo.afr = haplo[IID %in% afr.covar$IID]

# plink --bfile
# ~/mitonuclear2/data/genotype/MT/PMBB_Release-2020-2.0.${pop}.MTall \ --glm
# hide-covar \ --pheno
# ~/mitonuclear2/data/phenotype/PMBB_plink4mat_phecodes_.pheno \ --covar
# ~/mitonuclear2/data/phenotype/PMBB_plink4mat_AFR.covar \ --out
# ~/mitonuclear2/data/mtgwas/snp/AFR/pmbb_mtgwas_phecodes.AFR \ plink \ --bfile
# ~/mitonuclear2/data/genotype/MT/PMBB_Release-2020-2.0.EUR.MTall \ --glm
# hide-covar \ --pheno
# ~/mitonuclear2/data/phenotype/PMBB_plink4mat_phecodes_EUR.pheno \ --covar
# ~/mitonuclear2/data/phenotype/PMBB_plink4mat_EUR.covar \ --out
# ~/mitonuclear2/data/mtgwas/snp/EUR/pmbb_mtgwas_phecodes.EUR
