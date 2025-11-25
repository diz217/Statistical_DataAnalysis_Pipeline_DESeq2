library(GEOquery)

save_gse_meta = function(gse_id,outdir='~/R space/pipelinev0') 
{
  gse = getGEO(gse_id,GSEMatrix = TRUE)
  if (is.list(gse)) 
  {
    eset = gse[[1]]
  }else
  {
    eset = gse
  }
  meta = pData(eset)
  pathe = file.path(outdir,paste0(gse_id,"_series.txt"))
  write.table(meta,pathe,sep='\t',quote=FALSE)
  message("meta saved to",pathe)
}
args = commandArgs(trailingOnly = TRUE)
gse_id = args[1]
outdir = args[2]
save_gse_meta(gse_id,outdir)