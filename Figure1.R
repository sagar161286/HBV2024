# Following codes were used to generate the Figure 1 of the manuscript.
# Download the file hbv_hcv_github.rds from here: https://drive.google.com/file/d/1x5v8Vhiw1kOv-sSnfY-0iSd5NXYjawl2/view?usp=sharing

library(Seurat)
library(ggplot2)

hbv_hcv_github <- readRDS("hbv_hcv_github.rds")

# Figure 1A
DimPlot(hbv_hcv_github, reduction = "umap", group.by = "seurat_clusters_new", label.box = T, label = T)

# Figure 1B
fig1.genes <- read.table("Fig1B_heatmap_genes")
DotPlot(hcv.hbv.integrated, features =fig1.genes$V1,cols = c("BrBG"), cluster.idents = F,group.by = c("seurat_clusters_new")) + RotatedAxis()

