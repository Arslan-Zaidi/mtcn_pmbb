#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(data.table)
  library(optparse)
}))

option_list <- list(
  make_option(c("-p", "--pop"), type = "character", help = "path to pop file (must have 2 colums - IID and POP)"),
  make_option(c("-i", "--id", type = "character", help = "path to vcf id file (in the same order as vcf file)")),
  make_option(c("-o", "--out"), type = "character", help = "path to output file with extension")
)

args <- parse_args(OptionParser(option_list = option_list))

id = fread(args$id, header = FALSE)
pop = fread(args$pop)

id = fread("~/mitonuclear2/data/rfmix/classes/pmbb.afr.1kg.ceuryi.10.phased.ids", header = FALSE)
pop = fread("~/mitonuclear2/data/rfmix/1kg_ceu.yri.pop")

colnames(id) = c("IID")
colnames(pop) = c("IID","POP")

id = merge(id, pop, by = "IID", all.x = TRUE, sort = FALSE)
id[, sno := 1:.N]

id[, class := "0"]
id[POP == "YRI", class := "1"]
id[POP == "CEU", class := "2"]
id[is.na(POP), class := "0"]

id2 = copy(id)
id[, haplotype := 1]
id2[, haplotype := 2]

id3 = rbind(id, id2)
id3 = id3[order(sno, haplotype)]

write.table( t(id3$class), args$out, sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
