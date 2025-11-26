library(Matrix)
mtx <- readMM("~/R space/GSE253822_RAW/GSM8028380_Org_UT_D5_matrix.mtx.gz")
rownames(mtx)
dim(mtx)
genes = read.delim('~/R space/GSE253822_RAW/GSM8028380_Org_UT_D5_features.tsv.gz',header=FALSE)
cells = read.delim('~/R space/GSE253822_RAW/GSM8028380_Org_UT_D5_barcodes.tsv.gz',header=FALSE)
dim(genes)
dim(cells)
colnames(genes)
head(genes)
head(genes$V3)
colnames(cells)
head(cells)

mtx2 <- readMM("~/R space/GSE253822_RAW/GSM8028383_Org_AAV2_D5_matrix.mtx.gz")
mtx3 <- readMM("~/R space/GSE253822_RAW/GSM8028382_Org_Spk_D5_matrix.mtx.gz")
mtx4 <- readMM("~/R space/GSE253822_RAW/GSM8028381_Org_AAV9_D5_matrix.mtx.gz")

mtx22 <- readMM("~/R space/GSE253822_RAW/GSM8028379_Org_AAV2_48h_matrix.mtx.gz")
mtx33 <- readMM("~/R space/GSE253822_RAW/GSM8028378_Org_Spk_48h_matrix.mtx.gz")
mtx44 <- readMM("~/R space/GSE253822_RAW/GSM8028377_Org_AAV9_48h_matrix.mtx.gz")

mtx11 <- readMM("~/R space/GSE253822_RAW/GSM8028376_Org_UT_48h_matrix.mtx.gz")

rownames(mtx)=rownames(mtx11)=rownames(mtx2)=rownames(mtx3)=rownames(mtx4)=rownames(mtx22)=rownames(mtx33)=rownames(mtx44)=genes$V2
expr_bulk <- cbind(
  ut_d5 = Matrix::rowMeans(mtx),
  ut_d2 = Matrix::rowMeans(mtx11),
  aav2_d5 = Matrix::rowMeans(mtx2),
  aav2_d2 = Matrix::rowMeans(mtx22),
  spk_d5  = Matrix::rowMeans(mtx3),
  spk_d2  = Matrix::rowMeans(mtx33),
  aav9_d5 = Matrix::rowMeans(mtx4),
  aav9_d2 = Matrix::rowMeans(mtx44)
)

library(limma)
E = log2(expr_bulk+1)
grp <- factor(c("ut_d5","ut_d2","aav2_d5","aav2_d2","spk_d5","spk_d2","aav9_d5","aav9_d2"))
design = model.matrix(~0+grp)
colnames(design) = levels(grp)
fit <- lmFit(E, design)

fit2 <- eBayes(contrasts.fit(fit, makeContrasts(wow = ut_d5 - ut_d2, levels=design)))
res <- topTable(fit2, coef="wow", number=Inf)   
res$GeneID = rownames(res)
res = res[order(res$GeneID),]
write.csv(res, "~/R space/gse253822_utd5_utd2.csv")