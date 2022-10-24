

library(data.table)
F = rprojroot::is_rstudio_project$make_fix_file()

# load PCs
covar = fread(F("data/gcta/AFR/PMBB-Release-2020-2.0_genetic_genotype.AFR.maf0.01.pca.eigenvec"))
colnames(covar) = c("FID", "IID", paste("PC", seq(1, 20, 1), sep = ""))

# load other covariats
quants = fread(F("data/phenotype/PMBB_2.2_quant_08112022.pheno"))
# remove rlrmtCN (we know it's normal)
quants = quants[variable != "rlrmtCN"]
quants[, `:=`(scaled_value, scale(value)), by = "variable"]
quants[variable == "Height" & scaled_value < -5, `:=`(value, NA)]
quants[variable == "Height" & scaled_value < -5, `:=`(scaled_value, NA)]
# also remove negative values
quants[value < 0, `:=`(value, NA)]
quants[value < 0, `:=`(scaled_value, NA)]
# also remove values that are crazy high (e.g. > 7 SD)
quants = quants[scaled_value < 7 & scaled_value > -7]
quants = na.omit(quants)

ast = quants[variable == "AST"]
ast = merge(covar, ast[, .(PMBB_ID, scaled_value)], by.x = "IID", by.y = "PMBB_ID")


fwrite(ast, F("data/gcta/AFR/pmbb_AFR.pcs.ast_09022022.qcovar"), sep = "\t")


#### add DARC genotype

darc.a = fread(F("data/genotype/DARC/pmbb_rs281477.AFR.raw"), header = TRUE)
darc.ad = fread(F("data/genotype/DARC/pmbb_rs2814778_dom.AFR.raw"), header = TRUE)

darc.a = darc.a[, c(1, 2, 7)]
darc.ad = darc.ad[, c(1, 2, 7, 8)]

darc.a = merge(covar, darc.a, by = c("FID", "IID"))
darc.ad = merge(covar, darc.ad, by = c("FID", "IID"))

fwrite(darc.a, F("data/gcta/AFR/pmbb_AFR.pcs.darc.a_09022022.qcovar"), sep = "\t",
  na = "NA", quote = FALSE)
fwrite(darc.ad, F("data/gcta/AFR/pmbb_AFR.pcs.darc.ad_09022022.qcovar"), sep = "\t",
  na = "NA", quote = FALSE)


## add APOL1 genotype
apol1 = fread(F("data/genotype/APOL1/PMBB-Release-2020-2.0.2_rs73885319_A_G.AFR.raw"),
  header = TRUE)
apol1 = apol1[, c(1, 2, 7)]
apol1 = merge(covar, apol1, by = c("FID", "IID"))

fwrite(apol1, F("data/gcta/AFR/pmbb_AFR.pcs.apol1_09062022.qcovar"), sep = "\t",
  na = "NA", quote = FALSE)


#### local ancestry hits
lanc.hits = fread(F("data/rfmix/admixmap/pmbb.afr.admixmap.pruned.hits.raw"))
lanc.hits = lanc.hits[, c(1, 2, 7:27)]

lanc.hits = merge(covar, lanc.hits, by = c("FID", "IID"))
fwrite(lanc.hits, F("data/gcta/AFR/pmbb_AFR.pcs.lanchits_09062022.qcovar"), sep = "\t",
  na = "NA", quote = FALSE)


##### blood prs
blprs = fread(F("data/genotype/blgwas/prs/blgwas_pmbb_09082022.qc.sscore"))
colnames(blprs) = c("FID","IID","CT","DOSAGE","BASO","EOS","HCT","HGB",	"LYM",	"MCH",	"MCHC","MCV","MONO",	"MPV","NEU",	"PLT",	"RBC",	"RDW",	"WBC")
blprs$CT = blprs$DOSAGE = NULL

blprs = merge(covar, blprs, by = c("FID","IID"))
blprs.reduced = blprs[, .(FID, IID, BASO, EOS, LYM, MONO, NEU, PLT)]

fwrite(blprs.reduced, F("data/gcta/AFR/pmbb_AFR.pcs.blprs_09082022.qcovar"), sep = "\t",
  na = "NA", quote = FALSE)

fwrite(blprs, F("data/gcta/AFR/pmbb_AFR.pcs.blprs_09082022.extended.qcovar"), sep = "\t",
    na = "NA", quote = FALSE)
