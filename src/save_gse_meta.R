library(GEOquery)

save_gse_meta = function(gse_id) 
{
  outdir='~/R space/'
  gse = getGEO(gse_id,GSEMatrix = TRUE)
  if (is.list(gse)) 
  {
    for (i in seq_along(gse))
    {
    eset = gse[[i]]
    meta = pData(eset)
    pathe = file.path(outdir,paste0(gse_id,"_series_",i,".txt"))
    write.table(meta,pathe,sep='\t',quote=FALSE)
    message("meta ",i," saved to",pathe)
    }
  }else
  {
    eset = gse
    meta = pData(eset)
    pathe = file.path(outdir,paste0(gse_id,"_series.txt"))
    write.table(meta,pathe,sep='\t',quote=FALSE)
    message("meta saved to",pathe)
  }
}
args = commandArgs(trailingOnly = TRUE)
gse_id = args[1]
save_gse_meta(gse_id)