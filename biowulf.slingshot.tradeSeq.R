
.libPaths("~/R.packages/")

library(slingshot, lib.loc="~/R.packages/")
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(SingleCellExperiment)
library(scater)



# read in data ----

sce = readRDS("./sce.(no iPSC, PHA).Leiden.var-dev.genes.RDS")

# choose file name base ----

filename.base = "sce.(no iPSC, PHA).Leiden.var-dev.genes.slingshot.start8.end.3.10"



# slingshot ----

# run on SCE
# The output is a SingleCellExperiment object with slingshot results incorporated
# all of the results are stored in a PseudotimeOrdering object, which is added to the colData of the original object
# and can be accessed via colData(sce)$slingshot
sce = slingshot(sce, reducedDim = "UMAP",
                clusterLabels = colData(sce)$seurat_clusters, # tenx.seurat$seurat_clusters
                stretch=0, # do not extend curves beyond end-points
                start.clus=8, end.clus=c(3,10))

saveRDS(sce, file = paste0("./",filename.base,".RDS"))
print(paste0("saved: ",filename.base,".RDS"))



# parallel computing ----

BPPARAM = BiocParallel::bpparam()
BPPARAM$workers = 6



# tradeSeq ----

# it may be better to run the script in two phases
# optimal K is picked based on phase 1
# the two phases also have different resource requirements
# comment slingshot and phase 1 when ready for phase 2



# phase 1 ----

library(tradeSeq)

# evaluate K
icMat = evaluateK(sce, parallel=TRUE, BPPARAM=BPPARAM) # parallel=TRUE, BPPARAM=BPPARAM

saveRDS(icMat, file = paste0("./",filename.base,"icMat.RDS"))
print(paste0("saved: ",filename.base,"icMat.RDS"))



# phase 2 ----

library(tradeSeq)

# fit GAM model
sce = fitGAM(sce, nknots=8, parallel=TRUE, BPPARAM=BPPARAM) # parallel=TRUE, BPPARAM = BPPARAM

saveRDS(sce, file = paste0("./",filename.base,".fitGAM.RDS"))
print(paste0("saved: ",filename.base,"fitGAM.RDS"))


