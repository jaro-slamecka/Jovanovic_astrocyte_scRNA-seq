

# load libraries ----

library(dplyr)
library(Seurat)
library(ggplot2)
library(extrafont)
library(viridis)
library(Cairo)



# load 10X datasets ----

basedir = "../CELLRANGER/"
samples = c("hPSC-WA09-Astro-2", "hPSC-FCDI-Astro", "Glutamate-neurons")



# read 10X data in ----

tenx.data = lapply(samples,
                   function(x) {
                     tenx.data = Read10X(data.dir=paste0(basedir,x,"/outs/filtered_feature_bc_matrix/"))
                     return(tenx.data)
                   }
)

# rename the list with sample names for an easy iteration over them
names(tenx.data) = samples

# initialize the Seurat object with the raw (non-normalized data) ----
# all samples in one list
tenx.seurat.list = lapply(samples,
                          function(x) {
                            tenx.seurat = CreateSeuratObject(counts=tenx.data[[x]], project=x, min.cells=3, min.features=200)
                            return(tenx.seurat)
                          }
)

# merge Seurat objects ----
tenx.seurat = merge(x = tenx.seurat.list[[1]], y = tenx.seurat.list[2:length(tenx.seurat.list)],
                    project="Astro", add.cell.ids=samples)

# check number of cells from each original Seurat object
table(tenx.seurat$orig.ident)



# PRE-PROCESS ----

tenx.seurat = NormalizeData(tenx.seurat)
tenx.seurat = FindVariableFeatures(tenx.seurat)
tenx.seurat = ScaleData(tenx.seurat)
tenx.seurat = RunPCA(tenx.seurat)
tenx.seurat = RunUMAP(tenx.seurat, dims=1:20)
tenx.seurat = FindNeighbors(tenx.seurat)
tenx.seurat = FindClusters(tenx.seurat, algorithm=4, resolution=0.4)
tenx.seurat = RunTSNE(tenx.seurat)



# UMAP/tSNE ----

# UMAP
UMAPPlot(tenx.seurat, pt.size=1.2, label=TRUE, label.size=10) +
  update_geom_defaults("text", list(family="Noto Sans Cond", size=36)) +
  scale_color_discrete(c=100, l=58) +
  theme(text=element_text(size = 36, family="Noto Sans Cond"),
        axis.text=element_text(size=24),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.text=element_text(size=28, angle=45),
        legend.key.width=unit(0.1,"in"),
        legend.spacing.x=unit(0.01,"in"), legend.spacing.y=unit(0.2,"in"),
        legend.direction="horizontal", legend.position="top", legend.justification="center") +
  guides(color=guide_legend(label.position="top", nrow=1, override.aes=list(size=5)))
ggsave(filename="umap.clusters.pdf", device=cairo_pdf, width=8.2, height=8.6)

# UMAP - color by sample
plot.colors = c(viridis_pal(begin=0.2,end=0.7, option="D")(4))
UMAPPlot(tenx.seurat, group.by="orig.ident", pt.size=1.2) +
  scale_color_discrete(c=100, l=58) +
  update_geom_defaults("text", list(family="Noto Sans Cond", size=36)) +
  theme(text=element_text(size = 36, family="Noto Sans Cond"),
        axis.text=element_text(size=24),
        axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_blank(),
        legend.text=element_text(size=28, angle=0),
        legend.key.width=unit(0.1,"in"),
        legend.spacing.x=unit(0.01,"in"), legend.spacing.y=unit(0.2,"in"),
        legend.direction="vertical", legend.position="right", legend.justification="center") +
  guides(color=guide_legend(label.position="right", ncol=1, override.aes=list(size=5)))
ggsave(filename="umap.color.by.sample.pdf", device=cairo_pdf, width=12, height=8.4)
rm(plot.colors)



# FEATURE PLOTS ----

draw.feature.plots = function(seurat.object = tenx.seurat,
                              genes.to.plot = NULL,
                              unit.size = 5, # unit size of each plot ~ inches for pdf, multiplier factor for png
                              plot.ncol = 4,
                              plot.device = "pdf",
                              plot.file.name=NULL,
                              use.ggsave = TRUE) {
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(extrafont)
  library(viridis)
  library(Cairo)
  library(patchwork)
  # create the file name
  if (is.null(plot.file.name)) {plot.file.name=paste0("feature.plot.", paste(genes.to.plot[1:3], collapse = "."))}
  else {plot.file.name=paste0("feature.plot.", plot.file.name)}
  # make a list of individual plots
  feature.plots = FeaturePlot(seurat.object, features=genes.to.plot, pt.size=1.2, min.cutoff="q9", order=TRUE, combine=FALSE)
  # modify the theme individually
  feature.plots = lapply(X=feature.plots,
                         FUN=function(x) {x + theme(text=element_text(family="Noto Sans Cond"),
                                                    plot.title=element_text(size=32),
                                                    axis.title.x=element_blank(), axis.title.y=element_blank(),
                                                    legend.key.size=unit(0.3,"in"), legend.text=element_text(size=18))})
  # open graphical device
  if (plot.device=="png" && use.ggsave==FALSE) {
    CairoPNG(file=paste0(plot.file.name,".png"),
             width=unit.size*plot.ncol*100,
             height=(unit.size*length(genes.to.plot))/plot.ncol*100)}
  if (plot.device=="pdf" && use.ggsave==FALSE) {
    CairoPDF(file=paste0(plot.file.name,".pdf"),
             width=unit.size*plot.ncol*1,
             height=(unit.size*length(genes.to.plot))/plot.ncol*1)}
  # print is important for plotting from within a function
  # otherwise the plot is not sent to an open device
  print(wrap_plots(feature.plots, ncol=plot.ncol))
  if (use.ggsave==FALSE) {dev.off()}
  # save plots with ggsave
  if (plot.device=="png" && use.ggsave==TRUE) {
    ggsave(filename=paste0(plot.file.name,".png"), device=png,
           width=unit.size*plot.ncol*100, height=(unit.size*length(genes.to.plot))/plot.ncol*100)}
  if (plot.device=="pdf" && use.ggsave==TRUE) {
    ggsave(filename=paste0(plot.file.name,".pdf"), device=cairo_pdf,
           width=unit.size*plot.ncol*1, height=(unit.size*length(genes.to.plot))/plot.ncol*1)}
}


# choose features

gene.group = list()
# top markers per cluster
gene.group$general = all.markers %>% group_by(cluster) %>% top_n(n=3, wt=avg_log2FC)
# manual, split by categories/clusters:
gene.group$split = list()
gene.group$split$neuron = c("MAP2","DCX","SYT1","SYT4","EOMES","NEUROD1","NEUROD6","ROBO1","NSG2","INA")
gene.group$split$glia.neurons = c("FABP7","NFIA","HOPX","S100B","GFAP","MAP2","MAPT","NSG2","STMN2","SYT1")


# draw plots
draw.feature.plots(seurat.object=tenx.seurat, genes.to.plot=gene.group$split$neuron,
                   unit.size=5, plot.ncol=5, plot.device="pdf", plot.file.name="neurons")
draw.feature.plots(seurat.object=tenx.seurat, genes.to.plot=gene.group$split$glia.neurons,
                   unit.size=5, plot.ncol=5, plot.device="pdf", plot.file.name="glia.neurons")



# VIOLIN PLOTS ----

draw.violin.plots = function(seurat.object = tenx.seurat,
                             genes.to.plot = NULL,
                             unit.size = 5, # unit size of each plot ~ inches for pdf, multiplier factor for png
                             plot.ncol = 4,
                             plot.device = "pdf",
                             plot.file.name=NULL,
                             use.ggsave=TRUE) {
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(extrafont)
  library(viridis)
  library(Cairo)
  library(patchwork)
  # create the file name
  if (is.null(plot.file.name)) {plot.file.name=paste0("violin.plot.", paste(genes.to.plot[1:3], collapse = "."))}
  else {plot.file.name=paste0("violin.plot.", plot.file.name)}
  # make a list of individual plots
  violin.plots = VlnPlot(seurat.object, features=genes.to.plot, pt.size=1.2, combine=FALSE) #+ NoLegend()
  # modify the theme individually
  violin.plots = lapply(X=violin.plots,
                        FUN=function(x) {x + NoLegend() + # NoLegend() has to be here (it's a Seurat function)
                            scale_fill_discrete(c=100, l=58) + # violin fill color - to make it the same as UMAP point colors
                            theme(text=element_text(family="Noto Sans Cond"),
                                  plot.title=element_text(size=32),
                                  axis.text.x=element_text(size=28), axis.text.y=element_text(size=28),
                                  axis.title.x=element_blank(), axis.title.y=element_blank()
                            )})
  # open graphical device
  if (plot.device == "png" && use.ggsave==FALSE) {
    CairoPNG(file=paste0(plot.file.name,".png"),
             width=unit.size*plot.ncol*100,
             height=(unit.size*length(genes.to.plot))/plot.ncol*100)}
  if (plot.device == "pdf" && use.ggsave==FALSE) {
    CairoPDF(file=paste0(plot.file.name,".pdf"),
             width=unit.size*plot.ncol*1,
             height=(unit.size*length(genes.to.plot))/plot.ncol*1)}
  # print is important for plotting from within a function
  # otherwise the plot is not sent to an open device
  print(wrap_plots(violin.plots, ncol=plot.ncol))
  if (use.ggsave==FALSE) {dev.off()}
  # save plots with ggsave
  if (plot.device=="png" && use.ggsave==TRUE) {
    ggsave(filename=paste0(plot.file.name,".png"), device=png,
           width=unit.size*plot.ncol*100, height=(unit.size*length(genes.to.plot))/plot.ncol*100)}
  if (plot.device=="pdf" && use.ggsave==TRUE) {
    ggsave(filename=paste0(plot.file.name,".pdf"), device=cairo_pdf,
           width=unit.size*plot.ncol*1, height=(unit.size*length(genes.to.plot))/plot.ncol*1)}
}

draw.violin.plots(tenx.seurat,
                  genes.to.plot=c("NFIA","HOPX","S100B",
                                  "GFAP","AQP4","HEPACAM",
                                  "TBR1","SYT1","GRIA2",
                                  "ROBO1","NRXN1","INA",
                                  "STMN2","NCOR2"),
                  plot.ncol=7)



# DOT PLOTS ----

# plot top 30 markers per cluster

lapply(levels(Idents(tenx.seurat)), function(x) {
  CairoPNG(file=paste0("dot.plot.cluster.",x,".png"),width=1200,height=1200)
  # CairoPDF(file=paste0("dot.plot.cluster",x,".pdf"),width=12,height=12)
  sel.features = all.markers %>% filter(cluster==x) %>% top_n(n=30, wt=avg_log2FC)
  plot = DotPlot(tenx.seurat, features=sel.features$gene) +
    coord_flip() +
    scale_x_discrete(limits=rev) +
    theme(text=element_text(family="Noto Sans Cond", size=28),
          axis.text.x = element_text(size=28), axis.text.y=element_text(size=28))
  print(plot)
  dev.off()
})



# MARKER analysis ----

# find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers = FindAllMarkers(tenx.seurat, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# display top 2 for each cluster
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
# write into a file
library(data.table)
fwrite(all.markers, file="all.markers.txt", sep="\t")



# SLINGSHOT ----

# run on Biowulf or locally on Linux

library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(RColorBrewer)



# prepare data for slingshot analysis ----

# subset based on the most highly variable genes - to be combined with deviance-based filtering of features, below
# tradeSeq because is computationally intensive
length(rownames(tenx.seurat))
tenx.seurat = FindVariableFeatures(tenx.seurat, nfeatures = 4000)
var.features = VariableFeatures(tenx.seurat)

# deviance-based feature selection using package scry
# kstreet13's (Kelly Street) GitHub recommendation on feature selection
# she is one of the authors of tradeSeq
# "since variance is generally correlated with mean expression level, selecting for highly variable genes can end up being roughly the same as selecting highly expressed genes, which aren't always the most informative"
library(scry)
# in the current version of Seurat, the function as.SingleCellExperiment isn't working properly
# except for with DietSeurat for which dimreducs argument has to be declared in order for embeddings to be kept
sce = as.SingleCellExperiment(DietSeurat(tenx.seurat, dimreducs = c("umap","pca")))
sce = devianceFeatureSelection(sce, assay="counts", sorted=TRUE)
CairoPDF(file="features ranked by deviance.pdf", width=12, height=3)
plot(rowData(sce)$binomial_deviance, type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance")
abline(v=4000, lty=2, col="red")
dev.off()
# finally, select the top 4000 genes, they are sorted so the following can be used
dev.features = names(sce[1:4000, ])
# make a union of variable and deviant genes
length(unique(c(var.features,dev.features)))
sce2 = sce[unique(c(var.features,dev.features)), ]
# save the SCE object
saveRDS(sce2, file="sce.var-dev.genes.RDS")



# Biowulf run ----

# the RDS file is now ready to be copied to Biowulf
# the script biowulf.slingshot.tradeSeq.R contains instructions to run slighshot and tradeSeq on the SCE object
# slingshot can be run locally on a Linux/Mac but tradeSeq is computationally intensive
# the results are then loaded back here for local exploration



# load RDS file with the SCE object from Biowulf ----

# install the required package
devtools::install_github("LTLA/TrajectoryUtils")

# read in the sce object created on Biowulf
sce = readRDS(file="sce.var-dev.genes.slingshot.start8.end.3.10.RDS")



# UMAP plots with pseudotimes and trajectory lines ----

lapply(grep(colnames(colData(sce)), pattern="slingPseudotime_", value=TRUE), function(X) {
  require(grDevices)
  require(slingshot)
  require(BUSpaRse)
  require(tidyverse)
  require(tidymodels)
  require(Seurat)
  require(scales)
  require(viridis)
  require(Matrix)
  require(SingleCellExperiment)
  require(scater)
  require(RColorBrewer)
  CairoPDF(file=paste0("slingshot.trajectories.var-dev.start8.end.3.10.",X,".pdf"), width=10, height=10)
  colors = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plotcol = colors[cut(sce[[X]], breaks=100)]
  plotcol[is.na(plotcol)] = "grey" # assign grey color to NA values - values that do not belong to the lineage
  plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, col='black')
  dev.off()
})

# save the pseudotime colors
# 11 is the maximum number of colors
CairoPDF(file="pseudotime.colors.pdf", width=10, height=3)
display.brewer.pal(11, "Spectral")
dev.off()



# tradeSeq ----

library(tradeSeq)
library(patchwork)

# load the SCE object processed on Biowulf containing fitGAM results
sce = readRDS(file="sce.var-dev.genes.slingshot.start8.end.3.10.fitGAM.RDS")

# association test,
ATres = associationTest(sce)

head(ATres)
library(data.table)
fwrite(ATres, file="tradeSeq association test results.tsv", sep="\t", row.names=TRUE, col.names=TRUE)

# discovering progenitor marker genes
startRes = startVsEndTest(sce)
oStart = order(startRes$waldStat, decreasing=TRUE)
# make an order of significant genes
sigGeneStart = names(sce)[oStart[1:200]]
# filter out ribosomal and mitochondrial genes
sigGeneStart = sigGeneStart[!sigGeneStart %in% c(grep(sigGeneStart, pattern="^MT-", value=TRUE),
                                                 grep(sigGeneStart, pattern="^RP[SL]", value=TRUE))]
# automatically plot top genes
plot.list = list()
plot.list = lapply(1:36, function(x) {
  plotSmoothers(sce, as.matrix(assays(sce)$counts), gene=sigGeneStart[x]) +
    ggtitle(label=sigGeneStart[x]) +
    theme(text=element_text(size=24)) +
    guides(color=guide_legend(label.position="right", ncol=1, override.aes=list(size=5)))
})
CairoPNG(file="tradeSeq progenitor genes.png", width=3600, height=3000)
wrap_plots(plot.list, nrow=6)
dev.off()

# discovering differentiated marker genes
endRes = diffEndTest(sce)
oEnd = order(endRes$waldStat, decreasing=TRUE)
# make an order of significant genes
sigGeneEnd = names(sce)[oEnd[1:200]]
# filter out ribosomal and mitochondrial genes
sigGeneEnd = sigGeneEnd[!sigGeneEnd %in% c(grep(sigGeneEnd, pattern="^MT-", value=TRUE),
                                           grep(sigGeneEnd, pattern="^RP[SL]", value=TRUE))]
# automatically plot top genes
plot.list = list()
plot.list = lapply(1:36, function(x) {
  plotSmoothers(sce, as.matrix(assays(sce)$counts), gene=sigGeneEnd[x]) +
    ggtitle(label=sigGeneEnd[x]) +
    theme(text=element_text(size=24)) +
    guides(color=guide_legend(label.position="right", ncol=1, override.aes=list(size=5)))
})
CairoPNG(file="tradeSeq differentiated genes.png", width=3600, height=3000)
wrap_plots(plot.list, nrow=6)
dev.off()

# compare progenitor genes and genes of differentiated cells
venn = jaro.draw.venn(group1=sigGeneStart,
                      group2=sigGeneEnd,
                      group1.name="progenitors",
                      group2.name="differentiated")

# plot selected genes
plot.list = list()
plot.list = lapply(c("VIM","GFAP","GPM6B","PRDX6","CD63",
                     "CD24","RACK1","IGFBP2","FTL","FTH1",
                     "CLU","MAP2","GAP43","TMSB10","STMN2"), function(x) {
                       plotSmoothers(sce, as.matrix(assays(sce)$counts), gene=x) +
                         ggtitle(label=x) +
                         theme(text=element_text(size=24)) +
                         guides(color=guide_legend(label.position="right", ncol=1, override.aes=list(size=5)))
                     })
CairoPNG(file="tradeSeq selected genes.png", width=2800, height=1500)
CairoPDF(file="tradeSeq selected genes.pdf", width=34, height=18)
wrap_plots(plot.list, nrow=3)
dev.off()



# Seurat INTEGRATION with external data ----

# Nowakowski et al., 2017: Spatiotemporal gene expression trajectories reveal developmental hierarchies of the human cortex



# download and prepare data ----

library(data.table)

mat = fread("https://cells.ucsc.edu/cortex-dev/exprMatrix.tsv.gz")
meta = read.table("https://cells.ucsc.edu/cortex-dev/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
krieg = CreateSeuratObject(counts = mat, project = "cortex-dev", meta.data=meta)

unique(krieg$WGCNAcluster)



# pre-process ----

krieg = NormalizeData(krieg)
krieg = FindVariableFeatures(krieg)
krieg = ScaleData(krieg)
krieg = RunPCA(krieg)
krieg = RunUMAP(krieg, dims=1:20)



# WGCNA cluster contains a group of cells with missing cluster assignment
# these can cause problems when plotting UMAP
which(krieg@meta.data$WGCNAcluster == "")
krieg@meta.data[krieg@meta.data$WGCNAcluster == "", ]$WGCNAcluster = "none-missing"

# set WGCNAcluster as the active identity
Idents(krieg) = "WGCNAcluster"



# integration ----

features = SelectIntegrationFeatures(object.list = c(tenx.seurat,krieg))

# find anchors
anchors = FindIntegrationAnchors(object.list=c(tenx.seurat,krieg), dims=1:30, max.features=800) # default dims=1:30 and max.features=200

# integrate data, this creates a new assay "integrated" in Assays (use this instead of RNA)
seurat.int = IntegrateData(anchors, dims = 1:30)

DefaultAssay(seurat.int) = "integrated"

# Run regular Seurat processing on integrated data
seurat.int = ScaleData(seurat.int)
seurat.int = RunPCA(seurat.int)
seurat.int = RunUMAP(seurat.int, dims=1:30) # dims=1:30



# UMAP of integrated data ----

levels(Idents(seurat.int))

# assign colors
plot.colors.1 = viridis_pal(begin=0.2,end=0.7, option="D")(11) # 11 clusters in VJ data
plot.colors = c(plot.colors.1, rep("gray",48)) # 48 clusters in Nowakowski data
plot.colors[c(12,13,50:54,59)] = c("red","orange",rep("maroon3",5),"maroon3") # highlight selected Nowakowski clusters
rm(plot.colors.1)

# plot UMAP
UMAPPlot(seurat.int, pt.size=1.2, label=TRUE, label.size=10) +
  update_geom_defaults("text", list(family="Noto Sans Cond", size=36)) +
  scale_color_manual(values=plot.colors) +
  theme(text=element_text(size=24, family="Noto Sans Cond"),
        axis.text=element_text(size=24),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.text=element_text(size=20, angle=0),
        legend.key.width=unit(0.1,"in"),
        legend.spacing.x=unit(0.01,"in"), legend.spacing.y=unit(0.2,"in"),
        legend.direction="horizontal", legend.position="right", legend.justification="center") +
  guides(color=guide_legend(label.position="right", ncol=1, override.aes=list(size=5)))
ggsave(filename="umap.VJ-Nowakowski.pdf", device=cairo_pdf, width=16.8, height=16.8)


