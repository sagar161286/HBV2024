# Following codes were used to generate the Figure 1 of the manuscript.
# Download the file hbv_hcv_github.rds from here: https://drive.google.com/file/d/1fNUoG_vSV0YLZS2Q5-kWEJ3li10Pjjf-/view?usp=sharing

library(Seurat)
library(ggplot2)

hbv_hcv_github <- readRDS("hbv_hcv_github.rds")

############################ Figure 1A ###############################

DimPlot(hbv_hcv_github, reduction = "umap", group.by = "seurat_clusters_new", label.box = T, label = T)

############################ Figure 1B ###############################

fig1.genes <- read.table("Fig1B_heatmap_genes")
DotPlot(hcv.hbv.integrated, features =fig1.genes$V1,cols = c("BrBG"), cluster.idents = F,group.by = c("seurat_clusters_new")) + RotatedAxis()

############################ Figure 1C ###############################

library(dplyr)
library(ggpubr)

metadata <- as.data.frame(hbv_hcv_github@meta.data)

# Tex signature (see genes in the uploaded file - tex_signature.rds)
metadata <- metadata %>% mutate( seurat_clusters_new=factor(seurat_clusters_new,levels=c("4", "1", "3", "2", "0")) )

ggviolin(metadata, x = "seurat_clusters_new", y = "exh_up_utz", fill = "seurat_clusters_new",
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = list(c("0","1"),c("0","2"),c("0","3"),c("0","4"),c("1","2"),c("1","3"),c("1","4"),c("2","3"),c("2","4"),c("3","4")), label = "p.format")+ # Add significance levels
  stat_compare_means(label.y = 160) 

# Tml signature (see genes in the uploaded file - tml_signature.rds)
metadata <- metadata %>% mutate( seurat_clusters_new=factor(seurat_clusters_new,levels=c("4", "3", "1", "2", "0")) )

ggviolin(metadata, x = "seurat_clusters_new", y = "mem_up_utz", fill = "seurat_clusters_new",
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = list(c("0","1"),c("0","2"),c("0","3"),c("0","4"),c("1","2"),c("1","3"),c("1","4"),c("2","3"),c("2","4"),c("3","4")), label = "p.format")+ # Add significance levels
  stat_compare_means(label.y = 160) 


# Memory signature (see genes in the uploaded file - memory_signature.rds)
metadata <- metadata %>% mutate( seurat_clusters_new=factor(seurat_clusters_new,levels=c("2", "4", "3", "0", "1")) )
ggviolin(metadata, x = "seurat_clusters_new", y = "mem_up_utz_wo_year", fill = "seurat_clusters_new",
         add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = list(c("0","1"),c("0","2"),c("0","3"),c("0","4"),c("1","2"),c("1","3"),c("1","4"),c("2","3"),c("2","4"),c("3","4")), label = "p.format")+ # Add significance levels
  stat_compare_means(label.y = 100) 

############################ Figure 1D ###############################

library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)


genes <- c("TCF7--chr5","LEF1--chr4","SELL--chr1","CCR7--chr17","CXCR3--chrX","ZEB1--chr10","ID3--chr1","CD28--chr2","KLRG1--chr12","CX3CR1--chr3","PRDM1--chr6","TNF--chr6","TBX21--chr17","BHLHE40--chr3","RORC--chr1","RORA--chr15","ZEB2--chr2","ID2--chr2","GZMA--chr5","GZMB--chr14","GNLY--chr2","PRF1--chr10","IFNG--chr12","CD69--chr12","TNFRSF9--chr1","CD38--chr4","HLA-DRB1--chr6","HLA-DRA--chr6","PDCD1--chr2","ENTPD1--chr10","CTLA4--chr2","LAG3--chr12","TIGIT--chr3","TRIB1--chr8","BATF--chr14","TOX--chr8","TOX2--chr20")
scaled_exp <- DotPlot(hbv_hcv_github, features = genes,cols = c("RdYlBu"), dot.scale = 6, cluster.idents = F,group.by = c("seurat_clusters_new")) + RotatedAxis()
scaled_exp <- scaled_exp$data[,c(3,4,5)]
scaled_exp$id <- paste0("C", scaled_exp$id)
scaled_exp <- as.data.frame(scaled_exp %>% pivot_wider(names_from = features.plot, values_from = avg.exp.scaled))
rownames(scaled_exp) <- scaled_exp$id
scaled_exp$id <- NULL

scaled_exp <- as.data.frame(t(scaled_exp))
annotation_col = data.frame(Clusters =c( "C0","C1", "C2","C3","C4"))
rownames(annotation_col) = c("C0","C1", "C2","C3","C4")
annotation_row = data.frame(class = factor(rep(c("Memory", "Effector", "Cytotoxic","Activated","Exhausted"), c(8, 10, 5,5,9))))
rownames(annotation_row) = rownames(scaled_exp)

ann_colors = list(Clusters = c(C0 = clust.col[1], C1 = clust.col[2], C2 = clust.col[3], C3 = clust.col[4], C4 = clust.col[5]),class = c(Memory = "#7570B3", Effector = "#E7298A", Cytotoxic = "#66A61E",Activated="darkred",Exhausted = "red"))
hmcol = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
pheatmap(scaled_exp, cluster_rows = F, cluster_cols = F,annotation_col = annotation_col,  annotation_colors = ann_colors, show_colnames = F,annotation_row = annotation_row ,gaps_row = c(8,18,23,28,37), legend_labels = F, color = hmcol)


############################ Figure 1E ###############################

DimPlot(hbv_hcv_github, group.by = c("infection"))

set.seed(123)
clust.col = sample(rainbow(30))
types <- as.data.frame(hbv_hcv_github@meta.data)
types$new_types <- types$seurat_clusters_new
types2 <- types$new_types
counts <- as.data.frame.matrix(t(table(types2,types$infection)))
counts <- t(counts/rowSums(counts)*100)
rownames(counts) <- c(paste0("C",c(0:4)))
barplot(counts, col = clust.col, border = "white",ylim = c(0, 100), axes = FALSE)
axis(2, labels = c(0,20,40,60,80,100), at = c(0,20,40,60,80,100))

############################ Figure 1F ###############################

types <- as.data.frame(hbv_hcv_github@meta.data)
types$new_types <- types$seurat_clusters_new
types2 <- types$new_types
counts <- as.data.frame.matrix(t(table(types2,types$ml_ex)))
counts <- t(counts/rowSums(counts)*100)
rownames(counts) <- c(paste0("C",c(0:4)))
barplot(counts, col = clust.col, border = "white",ylim = c(0, 100), axes = FALSE)
axis(2, labels = c(0,20,40,60,80,100), at = c(0,20,40,60,80,100))


