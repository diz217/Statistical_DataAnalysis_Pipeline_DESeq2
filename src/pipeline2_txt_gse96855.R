# Pipeline: expression matrix -> limma -> stat
# Name: tpm_matrix: Transcripts Per Million, a type of expression matrix, which is count_matrix normalized.
## Step 1: download the tpm matrix by hand
d = read.delim("~/R space/GSE96855_AAVwt_RPKM/GSE96855_AAVwt_RPKM.txt", header = TRUE, row.names = 1,check.names = FALSE)
#d = read.delim("~/R space/GSE96855_Rep2_Cap2_RPKM/GSE96855_Rep2_Cap2_RPKM.txt", header = TRUE, row.names = 1,check.names = FALSE) # the other half of the dataset, run them separately
colnames(d)
cols <- grep("^WT_(048|094)_\\d+", colnames(d), value=TRUE)
#cols <- grep("^R2C2_(048|094)_\\d+", colnames(d), value=TRUE)
E <- log2(as.matrix(d[, cols]) + 1)


## Step 2: limma = log2+eBayes
# design: design is special in each gse case
grp <- factor(ifelse(grepl("_048_", cols), "h48", "h94")) # 2 replica
design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
# design # run this line to see what design looks like 
library(limma)
fit <- lmFit(E, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(h94vsh48 = h94 - h48, levels=design)))
res <- topTable(fit2, coef="h94vsh48", number=Inf)   

write.csv(res, "~/R space/WT_94vs48hr.csv")
#write.csv(res, "~/R space/R2C2_94vs48hr.csv")
