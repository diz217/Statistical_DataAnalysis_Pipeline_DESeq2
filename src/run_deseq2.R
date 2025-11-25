library(DESeq2)
run_deseq2 = function(gse_id,outdir='~/R space/pipelinev0/') 
{
 f_count = paste0(outdir,gse_id,"_count_gsm_matrix.tsv") 
 f_meta = paste0(outdir,gse_id,"_meta_gsm_condition.tsv") 
 
 #count
 outdir = paste0(outdir,"DESeq2/")
 dir.create(outdir,showWarnings = FALSE)
 
 counts = read.table(f_count,sep='\t',header=TRUE,check.names=FALSE,stringsAsFactors = FALSE)
 rownames(counts) = counts[[1]]
 counts[[1]] = NULL
 
 counts = as.matrix(counts)
 storage.mode(counts) = "integer"
 
 #meta
 meta = read.table(f_meta,sep='\t',header=TRUE,check.names=FALSE,stringsAsFactors = FALSE)
 sample_ids = meta$geo_accession
 titles = meta$title

 #prep
 counts = counts[,sample_ids]
 coldata = data.frame(row.names = sample_ids,condition = factor(titles))
 counts_f = counts[rowSums(counts >= 1) >= 2,]

 #run deseq2
 dds = DESeqDataSetFromMatrix(countData = counts_f,colData = coldata,design = ~condition)
 dds = DESeq(dds)
 unique_cond = levels(factor(titles))

 #pair comparisons
 for (i in 1:(length(unique_cond)-1)){
   for (j in (i+1):length(unique_cond)){
     ref = unique_cond[i]
     case = unique_cond[j]
     contrast_name <- paste0(case, "_vs_", ref)
     res = results(dds,contrast=c("condition",case,ref))
     res_df = as.data.frame(res)
     res_df$Gene_id = rownames(res_df)
     outpath = file.path(outdir,paste0(gse_id,"_deseq2_",contrast_name,".csv"))
     write.csv(res_df,outpath,row.names=FALSE)
   }
 }
}
args = commandArgs(trailingOnly = TRUE)
gse_id = args[1]
outdir = args[2]
run_deseq2(gse_id,outdir)
