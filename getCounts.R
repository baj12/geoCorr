library(GEOquery)
library(Biobase)
library(foreach)
library(doParallel)
library(snow)
cores=detectCores()
cl <- makeCluster(rep("localhost", 12), type = "SOCK")
registerDoParallel(cl)



gdsIds = read.csv(file = "HS-malat1.txt", header = F)

gId = gdsIds[1,1]

#for gId in gdsIds
foreach (gId %in% gdsIds) %dopar% {
  library(GEOquery)
  library(Biobase)
  gds <- getGEO(paste0("GDS", gId),GSEMatrix=TRUE)
  est = GDS2eSet(gds,do.log2 = F)
  est
  pData(est)
  ex = exprs(est)
  
  # for columns / samples
  samp = 1 
  foreach( samp = 1:ncol(ex)) %dopar% {
    exCol = ex[,samp]
    exCol = exCol[!is.na(exCol)]
    exCol = exCol[exCol>0]
    
    ratios = data.frame(denom = character(), ratio = numeric())
    # for genes
    ngene = 1
    ratios =  foreach (ngene = 1:length(exCol), .combine = rbind) %do% {
      gName = names(exCol[ngene])[[1]]
      # we only look at ratios where the gene name name in the zaehler/counter is bigger than in the denominator
      zaehler = exCol[names(exCol) > gName]
      ratio = data.frame(denom = gName, ratio = zaehler / exCol[gName])
      # ratios = rbind(ratios, ratio)
    }
    saveRDS(object = ratios, file = paste0("GDS", gId, ".", samp, ".RDS" ))
  }
}
stopCluster(cl)
