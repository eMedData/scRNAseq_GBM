print("Load R libraries")
library(dplyr)
library(data.table)
library(Seurat)
library(patchwork)
library(ggplot2)
library(utils)
library(tidyverse)

sample <- "COH056"
sample0 <- "COH084"
sample1 <- paste(sample0, sample, sep = "_")
platform <- "10x"
prefix0 <- paste(sample0, platform, sep = "_")
prefix1 <- paste(sample1, platform, sep = "_")


if (TRUE) {stop("The value is TRUE, so the script must end here")
  } else { #continue the script
print("Script did NOT end!")
  }

print("load counts")
setwd(dir = "/Volumes/T7/em_bild2/GBM/")
#skip if analyzing second sample from same batch
counts.df <- fread(file = paste(prefix0, "read_counts.txt", sep = "_"), fill = T)
counts.df[1:5, 1:5]
dim(counts.df)
gene.symbol <- function(counts.df) {
  counts.df <- subset(counts.df, counts.df$"Gene Symbol" !="")
  counts.df <- unique(counts.df, by = "Gene Symbol")
  counts.df <- as.data.frame(counts.df)
  rownames(counts.df) <- counts.df$"Gene Symbol"
  counts <- counts.df[,4:ncol(counts.df)]

}
counts <-  gene.symbol(counts.df)
counts[1:5, 1:5]
dim(counts)
rm(counts.df)
#saveRDS(counts, file = paste(prefix0, "gene.symbol_read_counts.rds", sep = "_"))
write.table(counts, file = paste(prefix0, "gene.symbol_read_counts.txt", sep="_"), 
            sep="\t", quote=F)
#
setwd(dir = "/Volumes/T7/em_bild2/GBM/")
counts <- read.table(file = paste0(prefix3, 
                                   "_gene.symbol_read_counts.txt"),
                     sep = "\t", header = T, quote = "")

#load counts from here if analyzing subsequent sample from same batch
counts <- readRDS(file = paste(prefix0, "gene.symbol_read_counts.rds", sep = "_"))
counts[1:5, 1:5]
dim(counts)


print("loading metadata into list and bind rows")
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/cell_metadata/list/"))
annot_list1 <- list.files(pattern = paste0("cell_metadata.", sample, "_[1-3].txt"))
annot_list1
annot1 <- bind_rows(lapply(annot_list1, fread))
dim(annot1)
head(annot1, 5)
tail(annot1, 5)
rm(annot_list1)
write.table(annot1, file = paste(prefix1, "cell_metadata_raw_combined.txt", sep = "_"), 
            sep = "\t", quote = F, row.names = F)
saveRDS(annot1, file = paste(prefix1, "cell_metadata_raw_combined.rds", sep = "_"))
print("done")

print("isolate individual patients counts")
counts1 <- counts[match(annot1$"Cell", colnames(counts))]
dim(counts1)
counts1[1:5, 1:5]

setwd(dir = paste0("D:/em_bild2/GBM/", prefix1))
write.table(counts1, file = paste(prefix1, "gene.symbol_counts.txt", sep = "_"), 
            sep = "\t", quote = F)
saveRDS(counts1, file = paste(prefix1, "gene.symbol_counts.rds", sep = "_"))
print("done!")

print("match metadata to corresponding counts!")
annot1 <- annot1[match(colnames(counts1), annot1$"Cell")]
rownames(annot1) <- annot1$"Cell"
print("done!")

#png(filename = paste(prefix1, "cell_number.png", sep = "_"), width = 1024, height = 768)
ggplot(as.data.frame(table(annot1$"Sample")), aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill = "#1f1d66") +
  theme_classic(base_size = 18) + xlab("Sample") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 1.5)) +
  labs(title = "CELLS PER SAMPLE") +
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()

print("create seurat objects for each patient!")
seu1 <- CreateSeuratObject(counts = counts1, project = sample, 
                           meta.data = annot1, min.cells = 5, min.features = 100)
seu1
dim(seu1)
head(seu1@meta.data, 5)
tail(seu1@meta.data, 5)
print("done!")

#.f = function() {
print("saving seurat objects!")
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat/"))
saveRDS(seu1, file = paste(prefix1, "seurat_vst_cc.raw.rds", sep = "_"))
print("done!")
#}

print("start seurat pipeline")
#Take sampes through seurat pipeline to extract metadata at end for singleR
#Standard Pre-processing Workflow

print("store percent.mt in seu obj")
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seu1[["percent.mt"]] <- PercentageFeatureSet(seu1, pattern = "^MT-")
head(seu1@meta.data, 5)
print("done")

print("visualize qc metrics before filtering cells!")

#png(filename = paste(prefix1, "vlnplot_qc.metrics_nCount_before_filting.png", sep="_"), width=1024, height=786)
VlnPlot(seu1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "Sample")
#dev.off()

#FeatureScatter() is typically used to visualize feature-feature relationships,but can be used for anything calculated by the object, i.e. columns in
#object metadata, PC scores etc.

#png(filename = paste(prefix1, "feat.scatt_qc.metrics_nCount.percent.mt_before_filting.png", sep="_"), width=1024, height=768)
FeatureScatter(seu1, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample")
#dev.off()
#png(filename = paste(prefix1, "feat.scatt_qc.metrics_nCount.nFeature_before_filting.png", sep="_"), width=1024, height=768)
FeatureScatter(seu1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample")
#dev.off()

if (TRUE) {stop("The value is TRUE, so the script must end here")
  } else { #continue the script
print("Script did NOT end!")
  }

print("filtering cells!")
seu1 <- subset(seu1, subset = nCount_RNA > 100 & nCount_RNA < 1.5e4 & nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 20)
print("done!")

print("visualize qc metrics after filtering cells!")
#png(filename = paste(prefix1, "vlnplot_qc.metrics_nCount_after_filting.png", sep="_"), width=1024, height=786)
VlnPlot(seu1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "Sample")
#dev.off()

#png(filename = paste(prefix1, "feat.scatt_qc.metrics_nCount.percent.mt_after_filting.png", sep="_"), width=1024, height=768)
FeatureScatter(seu1, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample")
#dev.off()
#png(filename = paste(prefix1, "feat.scatt_qc.metrics_nCount.nFeature_after_filting.png", sep="_"), width=1024, height=768)
FeatureScatter(seu1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample")
#dev.off()

print("mark cell cycle genes!")
#Get cell cycle genes
setwd(dir = "D:/em_bild2/")
cc.genes <- readLines("regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

seu1 <- CellCycleScoring(seu1, s.features = s.genes, g2m.features = g2m.genes)
head(seu1@meta.data, 5)
print("done!")

print("NormalizeData!")
seu1 <- NormalizeData(seu1, normalization.method = "LogNormalize", scale.factor = 10000)
seu1[["RNA"]]@data[1:5,1:5]
print("done!")

print("find top variable features, identify top 10!")
seu1 <- FindVariableFeatures(seu1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seu1), 10)
top10
plot1 <- VariableFeaturePlot(seu1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T, xnudge = 0, ynudge = 0)
#png(filename = paste(prefix1, "top.var.feat.png", sep = "_"), width = 1024, height = 768)
#LabelPoints(plot = plot1, points = top10, repel = T, xnudge = 0, ynudge = 0)
#dev.off()
plot1 + plot2
plot2
print("done")

print("scaling/regression!")
all.genes <- rownames(seu1)
seu1 <- ScaleData(seu1, features = all.genes, vars.to.regress = c("percent.mt", "nCount_RNA", "S", "G2M"))
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat/"))
saveRDS(seu1, file = paste(prefix1, "Seurat_vst_cc.scaled_first_run.rds", sep="_"))
print("done!")

print("run PCA!")
seu1 <- RunPCA(seu1, features = VariableFeatures(object = seu1))
print("done!")

print("visualize PCs!")
print(seu1[["pca"]], dims = 1:4, nfeatures = 5)
VizDimLoadings(seu1, dims = 1:2, reduction = "pca")
#png(filename = paste(prefix1, "dimplot_dimplot_pca.png", sep = "_"), width = 1024, height = 768)
DimPlot(seu1, reduction = "pca", group.by = "Sample")
#dev.off()
#png(filename = paste(prefix1, "dimheatmap_pca.png", sep="_"), width=1024, height=768)
DimHeatmap(seu1, dims = 1:15, cells = 500, balanced = TRUE)
#dev.off()
#png(filename = paste(prefix1, "elbowplot.png", sep="_"), width=1024, height=768)
ElbowPlot(seu1)
#dev.off()

#print("clustree analysis/saving")
#helps determine resolution for FindClusters()
#png(filename = paste(prefix1, "Seurat_clustree.png", sep="_"), width=1024, height=1024)
#clustree(seu1, prefix = "RNA_snn_res.")
#dev.off()
#print("done")

print("stop here to analyze graphs to determine PC for downstream analyses")

if (FALSE) {stop("The value is TRUE, so the script must end here")
  } else { #continue the script
print("Script did NOT end!")
  }

print("dimension reduction!")
seu1 <- FindNeighbors(seu1, dims = 1:11)
seu1 <- FindClusters(seu1, resolution = 0.5)
seu1 <- RunUMAP(seu1, dims = 1:11)
seu1 <- RunTSNE(seu1, dims = 1:11)
print("done!")

print("plot umap and tsne!")
PCs <- "11"
DimPlot(seu1, reduction = "umap", group.by = "Sample", label = F, pt.size = 0.5, 
        cols = c("darkslategray1","darkorchid", "deeppink")) +
  labs(title = paste(prefix1, "Sample", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))
#"darkslategray1","darkslategray4",

DimPlot(seu1, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T, pt.size = 0.5) +
  labs(title = paste(prefix1, "clusters", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))

  

FeaturePlot(seu1, 
            features = c("AGER", "HMGB1", "S100A9", "S100A8", "S100A4", "S100B", "S10012"), 
            pt.size = 0.5, ncol = 3) &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

VlnPlot(seu1, 
        features = c("AGER","HMGB1", "S100A9", "S100A8", "S100A4", "S10012", "S100B"),
        ncol = 3, group.by = "RNA_snn_res.0.5") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

FeaturePlot(seu1, 
            features = c("AGER", "IL1B", "TNF", "IL6", "GLO1", "IL10"), 
            pt.size = 0.5, ncol = 3) &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

VlnPlot(seu1, 
        features = c("AGER", "IL1B", "TNF", "IL6", "GLO1", "IL10"),
        ncol = 3, group.by = "Sample") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

############## Individual markers #############  
#Astrocyte Markers
FeaturePlot(seu1, features = c("GFAP", "ALDH1L1", "ATP13A4", "SOX2"), pt.size = 0.25)
FeaturePlot(seu1, features = c("CD44", "IDH1", "PROM1", "NES"), pt.size = 0.25)
FeaturePlot(seu1, features = c("SLC25A18", "SLC39A12", "ALDH1L1", "RFX4"),
            pt.size = 0.25)
#Oligodendrocytes
FeaturePlot(seu1, features = c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11"), 
            ncol = 3, pt.size=0.25)

#Neuron Markers
FeaturePlot(seu1, features = c("MAP2", "NeuN", "PDGFD", "RYR3"), pt.size = 0.25)
FeaturePlot(seu1, features = c("UCHL1", "ENO2", "SPTAN1", "SNAP25"),
            order = T, pt.size = 0.25)
FeaturePlot(seu1, features = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "ENO2", "SNAP25"),
            order = T, pt.size = 0.25, ncol = 3)

FeaturePlot(integrated.2, features = c("PTPRC", "TYROBP",
                                   "GFAP", "EGFR", 
                                   "MBP", "TF"), 
            ncol = 3, pt.size = 0.25)  &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

#Macrophages
FeaturePlot(int_obj2, features = c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"), 
            ncol = 3, order=TRUE, pt.size=0.25)

#neutrophiles
FeaturePlot(int_obj2, features = c("CYP1B1", "G0S2", "ITGAX", "CD163"),
            order=TRUE, pt.size=0.25)
#Monocytes
FeaturePlot(int_obj2, features = c("CD14", "LYZ", "S100A9", "FCGR3A"),
            order=TRUE, pt.size=0.25)

FeaturePlot(integrated.2, features = c("CCL24", "TYROBP", "CSF1R",
                                   "S100A9", "LYZ", "AIF1"),
            ncol = 3, pt.size = 0.25) &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

FeaturePlot(int_obj2, features = c( "CD3D", "CD3E",
                                    "GZMB","KLRF1"),
            ncol = 2, pt.size = 0.25) &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))


#B cells
FeaturePlot(seu1, features = c("IGHM", "CD19"),
            order = T, pt.size = 0.25)
#T-cells
FeaturePlot(seu1,
            features = c("CD2", "CD3D", "CD3E", "CD3G"),
            order=TRUE, pt.size=0.25)
#NK cells
FeaturePlot(seu1,
            features = c("GZMB","KLRF1", "SH2D1B", "TYROBP"),
            order=TRUE, pt.size=0.25)

#Normal vs Cancer
FeaturePlot(seu1, features = c("PTPRC", "MKI67", "PCNA"), order=TRUE, pt.size=0.25)

print("WRITE METADATA")

out <- data.table(Cell.ID = colnames(x = seu1))
out_meta <- cbind(out, seu1@meta.data)

#write.table(out_meta, paste(prefix1,"cell_metadata.UMAPcluster0.txt", sep="_"), sep = "\t", quote = F, col.names = T)
#saveRDS(out_meta, paste(prefix1, "cell_metadata.UMAPcluster0.txt", sep="_"))
#print("add one column with cluster")
#paste("Cluster", annot1$RNA_snn_res.0.5, sep = "")
#cluster <- as.data.frame(paste("Cluster", out_meta$RNA_snn_res.0.5, sep = ""))
#colnames(cluster) <- "RNA_snn_res.0.5"
#cluster$Cell.ID <- out_meta$Cell.ID
#out_meta <- merge(out_meta, cluster, by = "Cell.ID")

setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat"))
write.table(out_meta, paste(prefix1, "cell_metadata.UMAPcluster_first_run.txt", sep="_"), 
            sep = "\t", quote = F, col.names = T)
saveRDS(out_meta, paste(prefix1, "cell_metadata.UMAPcluster_first_run.rds", sep="_"))

  print("write counts")
counts_gene.symbol_filt <- GetAssayData(object = seu1, slot = "counts")
counts_gene.symbol_filt <- as.matrix(counts_gene.symbol_filt)
write.table(counts_gene.symbol_filt, file = paste(prefix1, "gene_symbol.filtered.counts_first.run.txt", sep="_"), sep = "\t", quote = F)
saveRDS(counts_gene.symbol_filt, file = paste(prefix1, "gene_symbol.filtered.counts_first.run.rds", sep="_"))
saveRDS(seu1, file = paste(prefix1, "Seurat_vst_cc.final.seu.obj_first_run.rds", sep="_"))

rm(plot1)
rm(plot2)
print("END OF SCRIPT!")
#save environment (if computationally possible)
