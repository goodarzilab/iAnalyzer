#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
if("bestNormalize" %in% rownames(installed.packages()) == FALSE) {install.packages("bestNormalize",repos = "http://cran.us.r-project.org")}
suppressMessages(library(bestNormalize))

args = commandArgs(trailingOnly=TRUE)
metadata <- args[[1]]
#metadata <- '/rasis/Projects/Tools/iAnalyzer/data/classic/metadata.txt'

meta <- read_tsv(metadata)

indir <- dirname(metadata)
setwd(indir)

files <- paste0(meta[1] %>% pull(),".cnt")
print(files)

## Design formula
func <- args[[2]]
#func <- "~sample.type"
ref <- args[[3]]
#ref <- "lo"

weight.method <- args[[4]]
#weight.method <- "SES" #or "n" or "SE"

outfile <- args[[5]]

meta <- as.data.frame(meta)
covariate <- strsplit(func,'~')[[1]][[2]]
meta[[covariate]] <- factor(meta[[covariate]])
meta[[covariate]] <- relevel(meta[[covariate]], ref=ref)

## Generate the count table
datalist <- lapply(files, function(x){read.table(file=x,header=FALSE,col.names=c("sgRNA", sub(".cnt$", "", x)), check.names=F)})
m <- Reduce(function(...) merge(..., by=1, all = TRUE), datalist)
rownames(m) <- m[,1]
m <- m[,-1]
m[is.na(m)]<-0
colnames(m) <- meta[1] %>% pull()
write.table(m,file=gsub(".txt","_sgRNA_raw_counts.txt", outfile), sep="\t", quote=FALSE, col.names=NA )

## Filter the count table
levels <- unique(meta[[covariate]])
set1 <- (meta[[covariate]] == levels[1])
set2 <- (meta[[covariate]] == levels[2])
tmp<-m
tmp[tmp>0] <- 1
f <- m[rowSums(tmp[,set1])>1 & rowSums(tmp[,set2])>1,]
print(sprintf("Count matrix was filtered down to %d rows (initially %d).", dim(f)[[1]], dim(m)[[1]]))

mc <- as.matrix(f+0.5)
nm <- log2(mc %*% diag(1e6/colSums(mc)))
colnames(nm) <- colnames(m)

## DESeq2 analysis
fit <- t(apply(nm,1,function(x){summary(lm(x~meta[[covariate]]))$coefficients[2,]}))
colnames(fit) <- c("logFC", "SE", "t.value","p.value")
write.table(fit, file=gsub(".txt","_sgRNAs.txt", outfile), sep="\t", quote=FALSE, col.names=NA )

###use non.targeting guides to transform data
non.target <- data.frame(fit[grepl("non-targeting", row.names(fit)),])
(BNobject <- bestNormalize(non.target$logFC))
target <- data.frame(fit[! grepl("non-targeting", row.names(fit)),])
target$z.score <- predict(BNobject, target$logFC)

#z.normalize based on empirical null
#target$z.score <- (target$logFC-mean(non.target$logFC))/sd(non.target$logFC)

df <- as_tibble(target, rownames='sgRNA')
if (weight.method == "n"){
  df <- df %>% mutate(w.z.score=z.score, w=1)
}else if (weight.method == "SE"){
  df <- df %>% mutate(w.z.score=z.score/SE, w=1/SE^2)
}else if (weight.method == "SES"){
  df <- df %>% mutate(w.z.score=z.score*abs(logFC)/SE, w=(logFC/SE)^2)
}else{
  df <- df %>% mutate(w.z.score=logFC/SE^2, w=1/SE^2)
}
df <- df %>% mutate(gene=gsub("_\\S+$","",df$sgRNA))
res <- df %>% group_by(gene) %>% summarise(rho=mean(logFC), combined.z=sum(w.z.score)/sqrt(sum(w)))
res <- res %>% mutate(pv = 2*pnorm(-abs(combined.z)))
write.table(res, outfile, quote = F, sep="\t", row.names = T)

res$t <- as.factor(abs(res$rho)>0.5 & res$pv<0.01)
g <- ggplot(data=as.data.frame(res), aes(x=rho, y=-log10(pv), colour=t), alpha=0.75, size=1) +
  geom_vline(aes(xintercept=0), colour="black",linetype="dashed", size=0.5)+
  #xlim(c(-3, 3)) + ylim(c(0, 6)) +
  xlab(expression(log[2]~'fold-change (rho)')) + ylab(expression(-log[10]~paste(italic('p'),'-value')))+
  geom_point(alpha=0.75, size=1) +
  #scale_y_log10()+
  scale_colour_manual(values=c("#11008f", "#FF9100", "FF0000"))+theme_bw(30) +
  theme(panel.background = element_rect(colour = "black",size=2), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=12))
ggsave(g, filename=gsub(".txt","_volcano.pdf", outfile))

non.target$z.score <- BNobject$x.t
non.target$pv <- 2*pnorm(-abs(non.target$z.score))
non.target$t <- as.factor(abs(non.target$logFC)>0.5 & non.target$pv<0.01)
fdr <-  sum(as.logical(non.target$t), na.rm=T)/sum(as.logical(res$t), na.rm=T) #0.088
sprintf("Empirical FDR = %.3f",fdr)
