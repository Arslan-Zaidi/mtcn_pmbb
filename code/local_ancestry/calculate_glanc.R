#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(data.table)
  library(optparse)
  #library(R.utils)
}))

option_list <- list(
  make_option(c("-c", "--chrom", type = "character", help = "chromosome number")),
  make_option(c("-k","--clean"), action = "store_true", default=FALSE, help = "To clean the Viterbi based on ForwardBackward file or not [default]")
)

args <- parse_args(OptionParser(option_list = option_list))

#specify name of forward-backward file
chrom = args$chrom

message("reading viterbi file")
vitname = paste("~/mitonuclear2/data/rfmix/output/pmbb.afr.1kg.ceuryi.", chrom, ".0.Viterbi.txt", sep = "")
vit = fread( vitname, header = F)

#sample_file
sample_file=fread("~/mitonuclear2/data/rfmix/classes/pmbb.afr.1kg.ceuryi.10.phased.ids", header = F)
colnames(sample_file) <- "IID"

#create a new matrix to populate "cleaned" local ancestry data
vit2<-as.matrix(vit)

#specify number of snps
nsnps<-nrow(vit)

#no. of individuals
ninds<-ncol(vit2)/2

#the sample file has individuals from reference/ancestral groups as well. remove these
sample_file<-sample_file[!grep('NA',sample_file$IID),]

if(args$clean){
  message("opening connection to fb file")
  fbname = paste("~/mitonuclear2/data/rfmix/output/pmbb.afr.1kg.ceuryi.", chrom, ".0.ForwardBackward.txt", sep = "")

  con <- gzfile(fbname, "r",compression=6)


  message("cleaning ancestry file")
  # go over every SNP and determine which ancestry (1,2, or 3) has the highest posterior probability.
  # match this with the results in the viterbi file
  # correct if doesn't match
  # if none of the three have probability of > 0.9, replace ancestry with NA
  pb<-txtProgressBar(min=0,max=nrow(vit),styl=3)
  for(i in 1:nsnps){
  line=readLines(con,n=1)
  if(length(line)==0){break}else{
    line<-as.numeric(unlist(strsplit(line,split=" ")))
    isplit<-split(line, ceiling(seq_along(line)/2))
    for(j in 1:length(isplit)){
      a<-isplit[[j]]
      max.index<-which(a==max(a))
      if(max.index!=vit2[i,j]){vit2[i,j]<-max.index}else{
        if(max(a)<0.9){vit2[i,j]<-NA}
      }
    }
  }
  setTxtProgressBar(pb,i)
  }


  close(con)

}


message("calculating global ancestry")
odd.index <- seq(1, ncol(vit2), 2)
even.index <- seq(2, ncol(vit2), 2)
glanc.odd <- apply( vit2[, odd.index], 2, function(x){
  a <- length( which(x==1))/length(x)
  b <- length( which(x==2))/length(x)
  miss <- length( which(is.na(x) == "TRUE")) / length(x)
  return( cbind(a,b,miss) )
})

glanc.even <- apply( vit2[,even.index], 2, function(x){
  a <- length( which(x == 1))/length(x)
  b <- length( which(x == 2))/length(x)
  miss <- length(which(is.na(x) == "TRUE"))/length(x)
  return(c(a,b,miss))
})

glanc <- (glanc.odd + glanc.even)/2
glanc <- as.data.table(t(glanc))
colnames(glanc) <- c("African","European","unk")
glanc$chr <- chrom
glanc$nsnps <- nsnps
glanc$IID <- sample_file$IID

if(args$clean){
  fwrite(as.data.table(vit2),paste("~/mitonuclear2/data/rfmix/output/pmbb.afr.1kg.ceuryi.", chrom, ".0.Cleaned.Viterbi.txt",sep = ""),
          sep = "\t", quote = F, showProgress = T, col.names = F, na = "NA")

  fwrite(glanc,
    paste("~/mitonuclear2/data/rfmix/glanc/pmbb.afr.1kg.ceuryi.", chrom, ".0.Cleaned.glanc",sep = ""),
    sep = "\t", quote = F, showProgress = T,col.names = T, na = "NA")
}else{
  fwrite(glanc,
    paste("~/mitonuclear2/data/rfmix/glanc/pmbb.afr.1kg.ceuryi.", chrom, ".0.glanc",sep = ""),
    sep = "\t", quote = F, showProgress = T,col.names = T, na = "NA")
}
