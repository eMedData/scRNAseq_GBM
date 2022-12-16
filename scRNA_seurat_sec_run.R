print("Load R libraries")
library(dplyr)
library(data.table)
library(Seurat)
library(patchwork)
library(ggplot2)
library(utils)
library(RColorBrewer)
library(Matrix)

detach("package:dplyr", unload = T)

if (FALSE) {stop("The value is TRUE, so the script must end here")
  } else { #continue the script
print("Script did NOT end!")
  }

sample <- "UR1113"
sample0 <- "COH056"
sample1 <- paste(sample0, sample, sep = "_")
platform <- "10x"
prefix0 <- paste(sample0, platform, sep = "_")
prefix1 <- paste(sample1, platform, sep = "_")

c1 <-

print("load raw seurat obj from first run")
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat"))
seu1 <- readRDS(file = paste(prefix1, 
                             "seurat_vst_cc.raw.rds", sep="_"))
dim(seu1@meta.data)
head(seu1@meta.data, 10)
print("done")

print("loading metadata with singleR cell annotations")
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/singleR"))
annot <- read.table(paste(prefix1, "cell_metadata.UMAPcluster.annot_singleR.txt", sep="_"), 
                    sep = "\t", quote = "", header = T)
dim(annot)
head(annot, 10)

print("match raw seurat obj to singleR celltypes annotations")
dim(seu1@meta.data)
seu1@meta.data <- seu1@meta.data %>%
  filter(rownames(seu1@meta.data) %in% annot$Cell.ID)
dim(seu1@meta.data)
dim(seu1)
seu1 <- seu1[,colnames(seu1) %in% annot$Cell.ID]
dim(seu1)
head(seu1@meta.data, 10)

print("match annot to seu1/inset celltypes in obj")
annot <- annot[match(seu1$Cell, annot$Cell.ID),]
head(annot, 10)
seu1[["HPCA_labels_main"]] <- annot$HPCA_labels_main
seu1[["ENCODE_labels_main"]] <- annot$ENCODE_labels_main
dim(seu1@meta.data)
head(seu1@meta.data, 10)

print("doublets")
setwd("D:/em_bild2/GBM/doublet_gbm/output")
fluid1 <- read.table(paste(sample, "Fluid1.scrublet_out.csv", sep = "-"),
                     sep = "\t", header = T, quote = "")
fliud1 <- sort(fluid1$scrublet_doublet_call1, decreasing = T)

fluid2 <- read.table(paste(sample, "Fluid2.scrublet_out.csv", sep = "-"),
                     sep = "\t", header = T, quote = "")
fliud2 <- sort(fluid2$scrublet_doublet_call1, decreasing = T)

tumor_cent <- read.table(paste(sample, "TumorCenter.scrublet_out.csv", 
                               sep = "-"),
                         sep = "\t", header = T, quote = "")
tumor_cent <- tumor_cent[order(tumor_cent$scrublet_doublet_call1, 
                               decreasing = T),]

tumor_edge <- read.table(paste(sample, "TumorEdge.scrublet_out.csv",
                               sep = "-"),
                         sep = "\t", header = T, quote = "")
tumor_edge <- tumor_edge[order(tumor_edge$scrublet_doublet_call1, 
                               decreasing = T),]

#doublet <- rbind(fluid1, fluid2, tumor_cent, tumor_edge)
#doublet <- rbind(fluid1,tumor_cent, tumor_edge)
doublet <- rbind(tumor_cent, tumor_edge)
doublet <- doublet[order(doublet$scrublet_doublet_call1,
                         decreasing = T),]
doublet_out <- subset(doublet, doublet$scrublet_doublet_call1 != "False")

dim(seu1@meta.data)
doublet_out$Cell.ID %in% seu1@meta.data$Cell
seu1@meta.data <- seu1@meta.data[!(seu1@meta.data$Cell %in% doublet_out$Cell.ID),]
doublet_out$Cell.ID %in% seu1@meta.data$Cell
dim(seu1@meta.data)

dim(seu1)
doublet_out$Cell.ID %in% colnames(seu1)
seu1 <- seu1[,!(colnames(seu1) %in% doublet_out$Cell.ID)]
doublet_out$Cell.ID %in% colnames(seu1)
dim(seu1)

print("remove cells >50")
sort(table(seu1@meta.data$HPCA_labels_main), decreasing = T)
dim(seu1@meta.data)
Unknown <- seu1@meta.data %>%
  group_by(HPCA_labels_main) %>%
  tally() %>%
  filter(n<50)

seu1@meta.data <- seu1@meta.data %>%
  mutate(HPCA_labels_main = ifelse(is.na(HPCA_labels_main), 
                                   "Unknown", 
                                   HPCA_labels_main)) %>%
  mutate(HPCA_labels_main = ifelse(HPCA_labels_main %in% 
                                     Unknown$HPCA_labels_main, 
                                   "Unknown", 
                                   HPCA_labels_main))
sort(table(seu1@meta.data$HPCA_labels_main), decreasing = T)
seu1 <- subset(seu1, subset = HPCA_labels_main != "Unknown")
seu1 <- subset(seu1, subset = HPCA_labels_main != "Tissue_stem_cells")
seu1 <- subset(seu1, subset = HPCA_labels_main != "Neurons")
sort(table(seu1@meta.data$HPCA_labels_main), decreasing = T)
dim(seu1@meta.data)
dim(seu1)

print("seurat workflow")
seu1[["percent.mt"]] <- PercentageFeatureSet(seu1, pattern = "^MT-")

print("mark cell cycle genes!")
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
plot2
print("done")

print("scaling/regression!")
all.genes <- rownames(seu1)
seu1 <- ScaleData(seu1, features = all.genes, vars.to.regress = c("percent.mt", "nCount_RNA", "S", "G2M"))
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat/"))
#saveRDS(seu1, file = paste(prefix1, "Seurat_vst_cc.scaled_sec_run.rds", sep="_"))
seu1 <- readRDS(file = paste(prefix1, "Seurat_vst_cc.scaled_sec_run.rds", sep="_"))
print("done!")
seu1
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
?FindClusters
print("dimension reduction!")
seu1 <- FindNeighbors(seu1, dims = 1:19)
seu1 <- FindClusters(seu1, resolution = 0.5)
seu1 <- RunUMAP(seu1, dims = 1:19)
seu1 <- RunTSNE(seu1, dims = 1:19)
print("done!")

pcs <- "15"
head(seu1@meta.data, 5)
DimPlot(seu1, reduction = "umap", group.by = "HPCA_labels_main", label = T, 
        repel = T, pt.size = 0.5, label.size = 3.5) +
  labs(title = paste(prefix1, "HPCA_labels_main", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 13, face = "bold"))
DimPlot(seu1, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T, pt.size = 0.5,
        repel = T) +
  labs(title = paste(prefix1, "Clusters", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))
DimPlot(seu1, reduction = "umap", label = T, pt.size = 0.5, group.by = "Sample",
        repel = T, cols = c("darkslategray1","darkorchid", "deeppink"))
#,"darkslategray4"
FeaturePlot(seu1,reduction = "umap", features = c("PTPRC", "TYROBP",
                                       "GFAP", "EGFR","SLC1A3", 
                                       "SOX4","MBP", "TF","CDK18"),
            pt.size = 0.25, ncol = 3) &
  theme(axis.title = element_text(size = 10, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 14, face = "bold"))
#seu1@meta.data$HPCA_labels_main[seu1@meta.data$HPCA_labels_main == "Neurons"] <- "Astrocyte"
#seu1@meta.data$HPCA_labels_main[seu1@meta.data$HPCA_labels_main == "Tissue_stem_cells"] <- "Astrocyte"

FeaturePlot(seu1, reduction = "umap", features = c("GAD1", "DLX6-AS1",
                               "MEG3", "GOLGA8B", 
                               "CACNA1A", "SYT1", "ENO2"),
            pt.size = 0.25, ncol = 3)
FeaturePlot(seu1, reduction = "umap", features = c("GAD1", "DLX6-AS1",
                                                   "MEG3", "GOLGA8B", 
                                                   "CACNA1A", "SYT1", "ENO2"),
            pt.size = 0.25, ncol = 3)

FeaturePlot(seu1, reduction = "umap", features = c("NDRG2", "GFAP","S100B"))
c10 <- FindMarkers(seu1, ident.1 = 7, min.pct = 0.25, test.use = "t", only.pos = T)
head(c10, 20)



print("assign cell clusters according to where they cluster as needed")
options(ggrepel.max.overlaps = Inf)
seu2 <- seu1
seu2$RNA_snn_res.0.5 <- as.character(seu2$RNA_snn_res.0.5)

seu2@meta.data <- within(seu2@meta.data, 
                         RNA_snn_res.0.5[RNA_snn_res.0.5 == "6" & HPCA_labels_main == "T_cells"] <- "8")
seu2@meta.data <- within(seu2@meta.data, 
                         RNA_snn_res.0.5[RNA_snn_res.0.5 == "8" & HPCA_labels_main == ""] <- "11")
seu2@meta.data <- within(seu2@meta.data, 
                         RNA_snn_res.0.5[RNA_snn_res.0.5 == "6" & HPCA_labels_main == "NK_cell"] <- "11")
seu2@meta.data <- within(seu2@meta.data, 
                         RNA_snn_res.0.5[RNA_snn_res.0.5 == "5" & HPCA_labels_main == "Neutrophils"] <- "12")
seu2@meta.data <- within(seu2@meta.data, 
                         RNA_snn_res.0.5[RNA_snn_res.0.5 == "4" & HPCA_labels_main == "Astrocyte"] <- "2")

seu2$RNA_snn_res.0.5 <- as.factor(seu2$RNA_snn_res.0.5)

c4 <- c("SYT1")


print("change Idents")
head(Idents(seu1), 5)
#seu2[["old.idents"]] <- Idents(object = seu2)
Idents(seu2) <- seu2@meta.data$RNA_snn_res.0.5
head(Idents(seu2),5)
#seu2[["old.idents_RNA_snn_res.0.5"]] <- Idents(object = seu2)
seu2[["old.idents"]] <- Idents(object = seu2)
seu2 <- RenameIdents(object = seu2,
                     '0' = "clone1",
                     '1' = "clone1",
                     '2' = "clone1",
                     '3' = "clone1",
                     '4' = "Monocytes",
                     '5' = "clone1",
                     '6' = "clone1",
                     '7' = "T-cells",
                     '8' = "Macrophages",
                     '9' = "clone1",
                     '10' = "NK-cells",
                     '11' = "Oligodendrocytes",
                     '12' = "clone2")

seu2[["cell_idents"]] <- Idents(object = seu2)
head(Idents(seu2), 5)

DimPlot(seu2, reduction = "umap", group.by = "cell_idents", label = F, 
        repel = T, pt.size = 0.5, label.size = 3.5) +
  labs(title = paste(prefix1, "Cell Idents", sep = " ")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 13, face = "bold"))

DimPlot(seu2, reduction = "umap", group.by = "Sample", label = F, pt.size = 0.5, 
        cols = c("darkslategray1","darkorchid", "deeppink")) +
  labs(title = paste(prefix1, "Sample", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))
# "darkslategray4",

DimPlot(seu2, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T, pt.size = 0.5,
        repel = T) +
  labs(title = paste(prefix1, "Clusters", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))


print("change celltypes according to where they cluster as needed")
#seu1@meta.data <- within(seu1@meta.data, 
#                        cell_idents[cell_idents == "NK-cells" & RNA_snn_res.0.5 == "13"] <- "T-cells")
#seu1@meta.data <- within(seu1@meta.data, 
#                         RNA_snn_res.0.5[RNA_snn_res.0.5 == "3" & HPCA_labels_main == "Neutrophils"] <- "14")
#seu1@meta.data <- within(seu1@meta.data, 
#                         RNA_snn_res.0.5[RNA_snn_res.0.5 == "7" & HPCA_labels_main == "NK_cell"] <- "15")

seu1 <- seu2
sort(table(seu1@meta.data$HPCA_labels_main), decreasing = T)

sort(table(seu1@meta.data$cell_idents), decreasing = T)


FeaturePlot(seu1, features = c("AGER", "S100A9", "S100A8", "S100A4"), pt.size = 0.5) &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))



############## Individual markers #############  
#Astrocyte Markers
FeaturePlot(seu1, features = c("GFAP", "ALDH1L1", "EGFR", "SOX2"), pt.size = 0.25)
FeaturePlot(seu1, features = c("CD44", "IDH1", "PROM1", "NES"), pt.size = 0.25)
FeaturePlot(seu1, features = c("SLC25A18", "SLC39A12", "ALDH1L1", "RFX4"),
            pt.size = 0.25)
#Oligodendrocytes
FeaturePlot(seu1, features = c("MBP", "TF", "PLP1", "MAG", "MOG", "PDGFRA"), pt.size=0.25)

#Neuron Markers
FeaturePlot(seu1, features = c("MAP2", "CDK4", "PDGFD", "RYR3"), pt.size = 0.25)
FeaturePlot(seu1, features = c("UCHL1", "ENO2", "SPTAN1", "SNAP25"),
            order = T, pt.size = 0.25)
FeaturePlot(seu1, features = c("NEFL", "GABRA1", "SYT1", "SLC12A5", "ENO2", "SNAP25"),
            order = T, pt.size = 0.25, ncol = 3)

#Macrophages
FeaturePlot(seu1, features = c("AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"), 
            ncol = 3, order=TRUE, pt.size=0.25)

#neutrophiles
FeaturePlot(seu1, features = c("CYP1B1", "G0S2", "ITGAX", "CD163"),
            order=TRUE, pt.size=0.25)
FeaturePlot(seu1, features = c("PERM", "FCGR3A", "FCGR2B", "FCGR2A"),
            order=TRUE, pt.size=0.25)
#Monocytes
FeaturePlot(seu1, features = c("CD14", "LYZ", "S100A9", "FCGR3A"),
            order=TRUE, pt.size=0.25)

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


#################
print("WRITE METADATA")
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat"))

out <- data.table(Cell.ID = colnames(x = seu1))
out_meta <- cbind(out, seu1@meta.data)

write.table(out_meta, paste(prefix1, "cell_metadata.UMAPcluster_sec_run2.txt", sep="_"),
            sep = "\t", quote = F, col.names = T, row.names = F)
saveRDS(out_meta, paste(prefix1, "cell_metadata.UMAPcluster_sec_run2.rds", sep="_"))

print("cell_type.txt file from annot")
head(seu1@meta.data, 5)
dim(seu1@meta.data)
celltype <- seu1@meta.data[,c(4,23)]
dim(celltype)
head(celltype,10)
saveRDS(celltype, paste(prefix1, "cell_metadata.UMAPcluster.singleR_annot_cell_types2.rds", sep="_"))
write.table(celltype, paste(prefix1, "cell_metadata.UMAPcluster.singleR_annot_cell_types.txt", sep="_"), 
            sep = "\t", quote = F, row.names = F, col.names = F)

print("write counts")
counts_gene.symbol_filt <- GetAssayData(object = seu1, slot = "counts")
counts_gene.symbol_filt <- as.matrix(counts_gene.symbol_filt)
write.table(counts_gene.symbol_filt, file = paste(prefix1, "gene_symbol.filtered.counts_sec_run.txt", sep="_"), sep = "\t", quote = F)
saveRDS(counts_gene.symbol_filt, file = paste(prefix1, "gene_symbol.filtered.counts_sec_run.rds", sep="_"))
saveRDS(seu1, file = paste(prefix1, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run.rds", sep="_"))

rm(seu2)
save.image(file = paste(prefix1, "sec_run_environment.RData", sep = "_"))
print("END OF SCRIPT!")


setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat"))
#saveRDS(seu1, file = paste(prefix1, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run2.rds", sep="_"))
seu1 <- readRDS(file = paste(prefix1, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run2.rds", sep="_"))

print("cell_type.txt file from annot")
head(seu1@meta.data, 5)
dim(seu1@meta.data)
celltype <- seu1@meta.data[,c(4,23)]
dim(celltype)
head(celltype,10)
saveRDS(celltype, paste(prefix1, "cell_metadata.UMAPcluster.singleR_annot_cell_types2.rds", sep="_"))
write.table(celltype, paste(prefix1, "cell_metadata.UMAPcluster.singleR_annot_cell_types2.txt", sep="_"), 
            sep = "\t", quote = F, row.names = F, col.names = F)

