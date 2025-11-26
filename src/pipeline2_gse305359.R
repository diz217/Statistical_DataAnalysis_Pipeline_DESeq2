# Pipeline: expression matrix -> limma -> stat
# Name: tpm_matrix: Transcripts Per Million, a type of expression matrix, which is count_matrix normalized.
## Step 1: download the tpm matrix by hand
d = read.csv("~/R space/GSE305359_tpm_matrix/GSE305359_tpm_matrix.csv", check.names=FALSE)
#colnames(d) # check the column names, geneID is called 'ENSEMBL'
rownames(d) <- make.unique(d$ENSEMBL)  # somehow geneID is not unique. so we make it unique 
#head(rownames(d)) # check the rows are indeed geneId
expr <- as.matrix(d[ , c("PM1","PM2","PM3","NT1","NT2","NT3")]) 
E = log2(expr+1)

## Step 2: limma = log2+eBayes
# design: design is special in each gse case
grp <- factor(c(rep("PM",3), rep("NT",3))) # 6 samples -> 2 levels: PM and NT
design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
# design # run this line to see what design looks like 
library(limma)
fit <- lmFit(E, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(PMvsNT = PM - NT, levels=design)))
res <- topTable(fit2, coef="PMvsNT", number=Inf)   
# colnames(res) # log2FoldChange、P.Value、adj.P.Val, avg,B(log odds)
write.csv(res, "~/R space/PM_NT.csv")
