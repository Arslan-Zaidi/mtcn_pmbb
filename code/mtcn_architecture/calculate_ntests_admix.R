#!/usr/bin/env Rscript

library(coda)

ninds = 9324
nswitches1 = matrix(NA, nrow = 22, ncol = 9324)
nswitches2 = matrix(NA, nrow = 22, ncol = 9324)
for (chrom in c(3:5, 8:22)) {
  vit = as.matrix(fread(paste("~/mitonuclear2/data/rfmix/output/pmbb.afr.1kg.ceuryi.",
    chrom, ".0.Cleaned.Viterbi.txt.gz", sep = "")))



  odd.ix = seq(1, nhaps, 2)
  even.ix = seq(2, nhaps, 2)

  print("converting haplotype ancestry to dosage")
  vit[vit == 1] = 0
  vit[vit == 2] = 1

  vit2 = vit[, odd.ix] + vit[, even.ix]
  vit2 = t(vit2)

  nswitches1[chrom, ] = apply(vit2, 1, function(x) {
    n1 = diff(x[!is.na(x)])
    return(sum(abs(n1)))
  })
  nswitches2[chrom, ] = apply(vit2, 1, function(x) {
    n1 = spectrum0.ar(x[!is.na(x)])$spec
    return(n1)
  })
  print(chrom)
}

mean(apply(nswitches1, 2, function(x) {
  sum(x, na.rm = TRUE)
}))
# 118.5

mean(apply(nswitches2, 2, function(x) {
  sum(x, na.rm = TRUE)
}))
# 17821.33
