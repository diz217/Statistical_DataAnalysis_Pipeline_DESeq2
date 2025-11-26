install.packages('Matrix')
library(Matrix)
mtx <- readMM("~/R space/GSE256350_RAW/GSM8094159_AAV9_1_matrix.mtx.gz")
genes = read.delim('~/R space/GSE256350_RAW/GSM8094159_AAV9_1_features.tsv.gz',header=FALSE)
cells = read.delim('~/R space/GSE256350_RAW/GSM8094159_AAV9_1_barcodes.tsv.gz',header=FALSE)

mtx2 <- readMM("~/R space/GSE256350_RAW/GSM8094160_AAV9_2_matrix.mtx.gz")
mtx3 <- readMM("~/R space/GSE256350_RAW/GSM8094161_incipient_1_matrix.mtx.gz")
mtx4 <- readMM("~/R space/GSE256350_RAW/GSM8094162_incipient_2_matrix.mtx.gz")

rownames(mtx)  = rownames(mtx2)  = rownames(mtx3)  = rownames(mtx4) = genes$V1

expr_bulk <- cbind(
  AAV9_1 = Matrix::rowMeans(mtx),
  AAV9_2 = Matrix::rowMeans(mtx2),
  INC_1  = Matrix::rowMeans(mtx3),
  INC_2  = Matrix::rowMeans(mtx4)
)
library(limma)
E = log2(expr_bulk+1)
grp <- factor(c("aav9","aav9","incipient","incipient"))
design = model.matrix(~0+grp)
colnames(design) = levels(grp)
fit <- lmFit(E, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(aav9vsinc = aav9 - incipient, levels=design)))
res <- topTable(fit2, coef="aav9vsinc", number=Inf)   
res$GeneID = rownames(res)
res = res[order(res$GeneID),]
write.csv(res, "~/R space/gse256350_aav9_inc.csv")
