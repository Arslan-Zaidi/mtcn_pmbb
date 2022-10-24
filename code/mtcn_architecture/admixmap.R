library(data.table)

args = commandArgs(TRUE)
chrom = args[1]

print("reading viterbi file")
vit = as.matrix(fread(paste("~/mitonuclear2/data/rfmix/output/pmbb.afr.1kg.ceuryi.",
  chrom, ".0.Cleaned.Viterbi.txt.gz", sep = "")))

nhaps = ncol(vit)

odd.ix = seq(1, nhaps, 2)
even.ix = seq(2, nhaps, 2)

print("converting haplotype ancestry to dosage")
vit[vit == 1] = 0
vit[vit == 2] = 1

vit2 = vit[, odd.ix] + vit[, even.ix]
vit2 = t(vit2)

sample.ids = fread("~/mitonuclear2/data/rfmix/pmbb.afr.samples", header = FALSE)
colnames(sample.ids) = "IID"

print("processing phenotype files")
covar = fread("~/mitonuclear2/data/genotype/pmbb.afr.weighted.04222022.glanc")
covar = covar[, .SD, .SDcols = c("IID", "African")]
pheno = fread("~/mitonuclear2/data/phenotype/PMBB_plink4mat_combined.blood_residuals_07272022.pheno")

pheno = merge(pheno, covar, by = "IID")
pheno$sno = seq(1, nrow(pheno))

geno = fread("~/mitonuclear2/data/genotype/DARC/pmbb_rs2814778_dom.AFR.raw2")
geno$IID = geno$FID
pheno = merge(pheno, geno, by = c("FID", "IID"))

sample.ids = merge(sample.ids, pheno, by = "IID", sort = FALSE, all.x = TRUE)

pc.form = paste("PC", seq(1, 20), sep = "")
pc.form = paste(pc.form, collapse = "+")

form1 = paste("rlrmtcn ~ African", sep = "")
form2 = paste("rlrmtcn ~ African + rs2814778_T+rs2814778_HET", sep = "")

print("running admixture mapping")
nsnps = ncol(vit2)
mat1 = matrix(NA, nrow = nsnps, ncol = 4)
mat2 = matrix(NA, nrow = nsnps, ncol = 4)
for (i in 1:nsnps) {
  l1 = lm(data = sample.ids, form = paste(form1, "+vit2[,", i, "]", sep = ""))
  l2 = lm(data = sample.ids, form = paste(form2, "+vit2[,", i, "]", sep = ""))

  mat1[i, ] = summary(l1)$coefficients[3, ]
  mat2[i, ] = summary(l2)$coefficients[5, ]
  # print(i)
}

snps = fread(paste("data/rfmix/snplist/pmbb.afr.1kg.ceuryi.", chrom, ".phased.snplist",
  sep = ""))
colnames(snps) = c("chrom", "bp", "cm")

mat1 = as.data.table(mat1)
mat2 = as.data.table(mat2)
colnames(mat1) = colnames(mat2) = c("beta", "se", "tstat", "pvalue")
mat1 = cbind(snps, mat1)
mat2 = cbind(snps, mat2)

print("printing output")
fwrite(mat1, paste("~/mitonuclear2/data/rfmix/admixmap/pmbb.afr.", chrom, ".mtcn.admixmap.txt",
  sep = ""), sep = "\t")

fwrite(mat2, paste("~/mitonuclear2/data/rfmix/admixmap/pmbb.afr.", chrom, ".mtcn.admixmap.rs2814778.txt",
  sep = ""), sep = "\t")
