library(limma)
run_limma_v2 = function(gse_id,f_count)
{
  indir ='~/R space/'
  f_meta = paste0(indir,gse_id,"_meta_gsm_condition.tsv") 
  f_pair = paste0(indir,gse_id,"_control_experiments.tsv") 
  
  #count
  outdir = paste0(indir,"Limma/")
  dir.create(outdir,showWarnings = FALSE)

  counts = read.table(f_count,sep='\t',header=TRUE,check.names=FALSE,stringsAsFactors = FALSE)
  rownames(counts) = counts[[1]]
  counts[[1]] = NULL
  
  counts = as.matrix(counts)
  storage.mode(counts) = "double"
  
  #meta
  meta = read.table(f_meta,sep='\t',header=TRUE,check.names=FALSE,stringsAsFactors = FALSE)
  sample_ids = meta$geo_accession
  titles = make.names(meta$condition)
  
  #prep
  coldata = data.frame(row.names = sample_ids,condition = factor(titles))
  counts_f = counts[rowSums(counts >= 1) >= 2,]
  
  #run limma 
  grp = factor(titles)
  design = model.matrix(~0+grp)
  colnames(design) = levels(grp)
  E = log(counts_f+1)
  fit = lmFit(E,design)
  
  #pair comparisons
  pcs = read.table(f_pair,sep='\t',header=FALSE,check.names=FALSE,stringsAsFactors=FALSE)
  colnames(pcs) = c("case","ref")
  pcs$case = make.names(pcs$case)
  pcs$ref = make.names(pcs$ref)
  for (k in seq_len(nrow(pcs))) {
    case = pcs$case[k]
    ref = pcs$ref[k]
    contrast_name <- paste0(case, "_vs_", ref)
    fit2 = contrasts.fit(fit,makeContrasts(contrasts = paste0(case,"-",ref),levels=design))
    res_df = data.frame(Gene_id = rownames(fit2),log2FC = as.numeric(fit2$coefficients),log2Mean = as.numeric(fit2$Amean))             
    outpath = file.path(outdir,paste0(gse_id,"_limma_",contrast_name,".csv"))
    write.csv(res_df,outpath,row.names=FALSE)
  }
}
args = commandArgs(trailingOnly = TRUE)
gse_id = args[1]
f_count = args[2]
run_limma_v2(gse_id,f_count)