

library(data.table)

# calculate mtDNA copy number from seq depth

dp = fread("~/mitonuclear2/data/dx_downloads/vcf/pmbb_1_45000.dp.gz")
sam = fread("~/mitonuclear2/data/dx_downloads/vcf/pmbb_1_45000.samples", header = FALSE)
colnames(dp) = c("position", sam$V1)

mdp = melt(dp, id.vars = c("position"))

mdp[value == ".", `:=`(value, NA)]
mdp$value = as.numeric(mdp$value)

#calculate average depth on mtDNA for each position
mdp.by.position = mdp[, .(mean.mt = mean(value, na.rm = TRUE), median = median(value, na.rm = TRUE)), by = "position"]

#write to file
fwrite(mdp.by.position, "~/mitonuclear2/data/phenotype/PMBB_mtdepth_byposition.txt", col.names = TRUE, sep = "\t")


mdp.summary = mdp[, .(mean.mt = mean(value, na.rm = TRUE)), by = "variable"]
colnames(mdp.summary) = c("IID", "mean.mt")

mdp2.summary = mdp[position < 2000 | position > 3500, .(mean.mt = mean(value, na.rm = TRUE)), by = "variable"]

colnames(mdp2.summary) = c("IID", "mean.mt")

idmap = fread("~/mitonuclear2/data/pmbb_45k_gwas_sample_mappings_post_qc_08042021_map.txt")
idmap = idmap[, c(1, 3)]

mdp.summary = merge(mdp.summary, idmap, by.x = "IID", by.y = "SampleName")
mdp2.summary = merge(mdp2.summary, idmap, by.x = "IID", by.y = "SampleName")

mdp.summary[, `:=`(FID, PMBB_ID)]
mdp2.summary[, `:=`(FID, PMBB_ID)]
mdp.summary$IID = NULL
mdp2.summary$IID = NULL

colnames(mdp.summary)[2] = colnames(mdp2.summary)[2] = "IID"


## load depth from autosomes

dp.auto = fread("~/mitonuclear2/data/genotype/PMBB_Release-2020-2.0-MT.sites.subsampled.from.exome.GL.DP")
dp.auto.sam = fread("~/mitonuclear2/data/genotype/PMBB_Release-2020-2.0-MT.sites.subsampled.from.exome.GL.samples", header = FALSE)
colnames(dp.auto) = c("chrom","position","ignore",dp.auto.sam$V1)
dp.auto$ignore = NULL

mdp.auto = melt(dp.auto, id.vars = c("chrom","position"))

mdp.auto[value == ".", `:=`(value, NA)]
mdp.auto$value = as.numeric(mdp.auto$value)

mdp.auto.summary = mdp.auto[, .(mean.auto = mean(value, na.rm = TRUE)), by = "variable"]
colnames(mdp.auto.summary) = c("IID", "mean.auto")
mdp.auto.summary[, FID := IID]


#### merge the two and calculate mtCN
dp.both = merge(mdp.summary, mdp.auto.summary, by = c("FID","IID"))
dp2.both = merge(mdp2.summary, mdp.auto.summary, by = c("FID","IID"))

dp.both[, mtcn := mean.mt/mean.auto]
dp2.both[, mtcn := mean.mt/mean.auto]

dp.both[, log.mtcn := log(mtcn)]
dp2.both[, log.mtcn := log(mtcn)]

dp.both[, mean.mtcn := (mtcn - mean(mtcn))/sd(mtcn)]
dp2.both[, mean.mtcn := (mtcn - mean(mtcn))/sd(mtcn)]

dp.both[, log.mtcn := (log.mtcn - mean(log.mtcn))/sd(log.mtcn)]
dp2.both[, log.mtcn := (log.mtcn - mean(log.mtcn))/sd(log.mtcn)]


#write to file
fwrite(dp.both[, .(FID, IID, mean.mtcn, log.mtcn)], "~/mitonuclear2/data/phenotype/PMBB_mtCN_unfiltered_45k_10212022.txt",
  sep = "\t")

fwrite(dp2.both[, .(FID, IID, mean.mtcn, log.mtcn)], "~/mitonuclear2/data/phenotype/PMBB_mtCN_filtered_45k_10212022.txt",
  sep = "\t")
