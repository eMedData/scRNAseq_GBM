library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(Seurat)
library(utils)
library(circlize)
library(fastcluster)
library(RColorBrewer)
#library(corrplot)
library(tidyverse)

#memory.size()
#memory.limit(30000)

#### Load 10x Seurat objects 

sample0 <- "COH056"
sample1 <- "COH084"
sample2 <- "COH056_COH084_integrated"
sample3 <- "COH056_NP4396"
sample4 <- "COH056_UR1113"
sample5 <- "COH084_JS4587"
sample6 <- "COH084_RY8558"
sample7 <- "COH084_YY2224"
platform <- "10x"

prefix0 <- paste(sample0, platform, sep="_")
prefix1 <- paste(sample1, platform, sep="_")
prefix2 <- paste(sample2, platform, sep="_")
prefix3 <- paste(sample3, platform, sep="_")
prefix4 <- paste(sample4, platform, sep="_")
prefix5 <- paste(sample5, platform, sep="_")
prefix6 <- paste(sample6, platform, sep="_")
prefix7 <- paste(sample7, platform, sep="_")


print("Load seurat objects")
setwd(dir = "D:/em_bild2/GBM/cca/seurat_objects")
int_obj1 <- readRDS(file = paste(prefix3, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run_infercnv.rds", sep = "_"))
int_obj2 <- readRDS(file = paste(prefix4, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run_infercnv.rds", sep = "_"))
int_obj3 <- readRDS(file = paste(prefix5, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run_infercnv.rds", sep = "_"))
int_obj4 <- readRDS(file = paste(prefix6, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run_infercnv.rds", sep = "_"))
int_obj5 <- readRDS(file = paste(prefix7, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run_infercnv.rds", sep = "_"))

setwd(dir = "D:/em_bild2/GBM/cca/seurat_objects2")
int_obj1 <- readRDS(file = paste(prefix3, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run2.rds", sep = "_"))
int_obj2 <- readRDS(file = paste(prefix4, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run2.rds", sep = "_"))
int_obj3 <- readRDS(file = paste(prefix5, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run2.rds", sep = "_"))
int_obj4 <- readRDS(file = paste(prefix6, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run2.rds", sep = "_"))
int_obj5 <- readRDS(file = paste(prefix7, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run2.rds", sep = "_"))

setwd(dir = "/Volumes/T7/em_bild2/GBM/cca/seurat_objects3")
int_obj1 <- readRDS(file = paste(prefix3, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run3_infercnv.rds", sep = "_"))
int_obj2 <- readRDS(file = paste(prefix4, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run3_infercnv.rds", sep = "_"))
int_obj3 <- readRDS(file = paste(prefix5, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run3_infercnv.rds", sep = "_"))
int_obj4 <- readRDS(file = paste(prefix6, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run3_infercnv.rds", sep = "_"))
int_obj5 <- readRDS(file = paste(prefix7, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run3_infercnv.rds", sep = "_"))



int_obj1
dim(int_obj1@meta.data)
head(int_obj1@meta.data, 5)
int_obj2
dim(int_obj2@meta.data)
head(int_obj2@meta.data, 5)
int_obj3
dim(int_obj3@meta.data)
head(int_obj3@meta.data, 5)
int_obj4
dim(int_obj4@meta.data)
head(int_obj4@meta.data, 5)
int_obj5
dim(int_obj5@meta.data)
head(int_obj5@meta.data, 5)

table(int_obj1@meta.data$cell_idents3)

table(int_obj2@meta.data$cell_idents3)

table(int_obj3@meta.data$cell_idents3)

table(int_obj4@meta.data$cell_idents3)

table(int_obj5@meta.data$cell_idents3)

#Idents(int_obj1) <- int_obj1@meta.data$RNA_snn_res.0.5
#Idents(int_obj2) <- int_obj2@meta.data$RNA_snn_res.0.5
#Idents(int_obj3) <- int_obj3@meta.data$RNA_snn_res.0.5
#Idents(int_obj4) <- int_obj4@meta.data$RNA_snn_res.0.5
#Idents(int_obj5) <- int_obj5@meta.data$RNA_snn_res.0.5

#rm(int_obj5)
setwd(dir = "/Volumes/T7/em_bild2/GBM/cca/seurat_objects3/")
integrated.2 <- readRDS(file = paste(prefix2,"cca_seurat_obj.rds", sep = "_"))
integrated.2
dim(integrated.2)

integrated.2 <- seu1
#saveRDS(integrated.2, file = paste(prefix2,"cca_seurat_obj.rds", sep = "_"))

#saveRDS(anchors.2, file = paste(prefix2,"cca_anchors.rds",sep = "_"))


anchors.2 <- FindIntegrationAnchors(object.list = list(int_obj1, int_obj2, int_obj3, 
                                                       int_obj4, int_obj5), dims = 1:30)
integrated.2 <- IntegrateData(anchorset = anchors.2, dims = 1:30)
print("done")
# Run the standard workflow for visualization and clustering
integrated.2 <- ScaleData(object = integrated.2, verbose = FALSE, assay = "RNA")
integrated.2 <- FindVariableFeatures(object = integrated.2, verbose = FALSE, assay = "RNA")
integrated.2 <- RunPCA(object = integrated.2, npcs = 30, verbose = FALSE, assay = "RNA")
ElbowPlot(integrated.2)
integrated.2 <- RunUMAP(object = integrated.2, reduction = "pca", dims = 1:11, assay = "RNA")
integrated.2 <- RunTSNE(object = integrated.2, reduction = "pca", dims = 1:11, assay = "RNA")
DimPlot(object = integrated.2, reduction = "umap", group.by = "HPCA_labels_main")
DimPlot(object = integrated.2, reduction = "umap", group.by = "cell_idents3", label = T)
DimPlot(object = integrated.2, reduction = "umap", group.by = "Sample")
#DimPlot(object = integrated.2, reduction = "umap")
PCs.rna <- "11"

# Run the standard workflow for visualization and clustering
DefaultAssay(object = integrated.2) <- "integrated"
integrated.2 <- ScaleData(object = integrated.2, verbose = FALSE)
integrated.2 <- RunPCA(object = integrated.2, npcs = 30, verbose = FALSE)
ElbowPlot(integrated.2)
integrated.2 <- RunUMAP(object = integrated.2, reduction = "pca", dims = 1:19)
integrated.2 <- RunTSNE(object = integrated.2, reduction = "pca", dims = 1:19)
PCs.integrated <- "19"



integrated.2 <- RenameIdents(object = integrated.2,
                     'Tumor-Astrocytes' = "Malignant Astrocytes",
                     'Tumor-Oligodendrocytes' = "Oligodendrocytes",
                     'Neurons' = "Neurons",
                     'Macrophages' = "Macrophages",
                     'Monocytes' = "Monocytes",
                     'Neutrophils' = "Neutrophils",
                     'T-cells' = "T-cells",
                     'NK-cells' = "NK-cells")

integrated.2[["cell_idents3"]] <- Idents(object = integrated.2)
head(Idents(integrated.2), 5)

DimPlot(object = integrated.2, reduction = "umap", group.by = "HPCA_labels_main",
        cols = c("firebrick1", "cadetblue1", "blue", "cadetblue4", "aquamarine","chartreuse4")) +
  #labs(title = paste(prefix2, "HPCA Labels")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))

####################
setwd(dir = "/Volumes/T7/em_bild2/GBM/cca/figures/")
pdf(file = paste(prefix2, "umap_cell.idents3.pdf",sep = "_"), 
    width = 5, height = 5)
DimPlot(object = integrated.2, reduction = "umap",
        label = F, 
        repel = T, 
        pt.size = .25,
        cols = c("firebrick1", "hotpink", "purple", "cyan",
                 "blue","green", "cadetblue4", "chartreuse4"))  +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold")) + 
  NoLegend()
dev.off()

pdf(file = paste(prefix2, "umap_sample.pdf",sep = "_"), 
    width = 5, height = 5)
DimPlot(object = integrated.2, reduction = "umap", group.by = "Sample",
        cols = c("darkslategray1", "darkorchid", "deeppink", 
                 "darkslategray1", "darkslategray4", "darkorchid", "deeppink", 
                 "darkslategray1", "darkorchid", "deeppink",
                 "darkslategray1", "darkslategray4", "darkorchid", "deeppink", 
                 "darkslategray1", "darkorchid", "deeppink")) +
  labs(title = "") +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) + 
  NoLegend()
dev.off()

pdf(file = paste(prefix2, "umap_sample.patients.pdf",sep = "_"), 
    width = 5, height = 5)
DimPlot(object = integrated.2, reduction = "umap", group.by = "Sample",
        cols = c("orange","orange","orange","pink","pink","pink","pink",
                 "darkolivegreen3","darkolivegreen3","darkolivegreen3","cyan","cyan","cyan","cyan",
                 "mediumorchid1","mediumorchid1","mediumorchid1")) +
  labs(title = "") +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  NoLegend()
dev.off()


#+
#  labs(title = paste(prefix2, "HPCA Idents"))

##For paper
pdf(filename = paste(prefix2, "umap_cell.idents3.pdf",sep = "_"), 
     width = 5.67, height = 5.04)
DimPlot(object = integrated.2, reduction = "umap", label = T, repel = T, pt.size = .5,
        cols = c("firebrick1", "hotpink", "purple", 
                 "cyan", "blue","green", "cadetblue4", "chartreuse4"
        ))  +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  NoLegend()
dev.off()
DimPlot(object = integrated.2, reduction = "umap", group.by = "Sample") +
  labs(title = paste(prefix2, "Sample", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))

DimPlot(object = integrated.2, reduction = "umap", group.by = "Sample", 
        cols = c("orange","orange","orange","pink","pink","pink","pink",
                 "blue","blue","blue","cyan","cyan","cyan","cyan",
                 "purple","purple","purple")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))

#######
DimPlot(object = integrated.2, reduction = "umap", group.by = "Sample",
        cols = c("orange","orange","orange","pink","pink","pink","pink",
                 "darkolivegreen3","darkolivegreen3","darkolivegreen3","cyan","cyan","cyan","cyan",
                 "mediumorchid1","mediumorchid1","mediumorchid1")) +
  labs(title = "") +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  NoLegend()
  

DimPlot(object = integrated.2, reduction = "umap", group.by = "Sample",
        cols = c("darkslategray1", "darkorchid", "deeppink", 
                 "darkslategray1", "darkslategray4", "darkorchid", "deeppink", 
                 "darkslategray1", "darkorchid", "deeppink",
                 "darkslategray1", "darkslategray4", "darkorchid", "deeppink", 
                 "darkslategray1", "darkorchid", "deeppink")) +
  labs(title = paste(prefix2, "Sample", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))

DimPlot(object = integrated.2, reduction = "umap", group.by = "Sample",
        cols = c("darkslategray1", "darkorchid", "deeppink", 
                 "darkslategray1", "darkslategray4", "darkorchid", "deeppink", 
                 "darkslategray1", "darkorchid", "deeppink",
                 "darkslategray1", "darkslategray4", "darkorchid", "deeppink", 
                 "darkslategray1", "darkorchid", "deeppink")) +
  labs(title = "") +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) + 
  NoLegend()

######
DimPlot(object = integrated.2, reduction = "umap", group.by = "RNA_snn_res.0.5",
        label = T) +
  labs(title = paste(prefix2, "Clusters", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))


DimPlot(object = integrated.2, reduction = "umap", group.by = "cell_idents3", 
        label = T, pt.size = 0.5) +
    labs(title = paste(prefix2, "Seurat Clusters", sep = "_")) +
    theme(axis.title = element_text(size = 10, face = "bold")) +
    theme(axis.line = element_line(size = 1.5)) +
    theme(axis.text = element_text(size=10, face = "bold")) +
    theme(legend.text = element_text(size = 14, face = "bold"))


DefaultAssay(object = integrated.2) <- "RNA"


#GBM Markers
#AC-LIKE

FeaturePlot(integrated.2, features = c("PTPRC", "TYROBP",
                                       "GFAP", "EGFR", "SLC1A3", 
                                       "ENO2", "TF", "MBP"),
                                       pt.size = 0.25, ncol = 3) &
  theme(axis.title = element_text(size = 10, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 14, face = "bold"))
  
FeaturePlot(integrated.2, features = c("GFAP", "EGFR", "ALDH1L1",
                                       "TF", "MBP", "PLP1", 
                                       "SYT1", "MAP2", "ENO2"),
            pt.size = 0.25, ncol = 3) &
  theme(axis.title = element_text(size = 10, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 14, face = "bold"))

FeaturePlot(integrated.2, features = c("GFAP", "CST3", "S100B", 
                               "EGFR", "SLC1A3", "HOPX"), 
            ncol = 3, pt.size = 0.25)
#MES-LIKE
FeaturePlot(integrated.2, features = c("DDIT3", "ENO2", "VIM", "HILPDA"), pt.size = 0.25)

#OLIGO-LIKE
FeaturePlot(integrated.2, features = c("MBP", "TF", "PLP1", "MAG", "MOG", "CLDN11"), 
            ncol = 3, pt.size=0.25)
#NPC-LIKE
FeaturePlot(integrated.2, features = c("STMN4", "STMN2", "DLX5", "SOX11"),
            pt.size = 0.25)


#immune
FeaturePlot(integrated.2, features = c("PTPRC", "TYROBP", "ITGAX", "CD163",
                                       "FCGR3A","LYZ", "S100A9", 
                                       "CD3D","GZMB"),
            pt.size = 0.25, ncol = 3) &
  theme(axis.title = element_text(size = 10, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 14, face = "bold"))


head(integrated.2@meta.data, 5)


DimPlot(integrated.2, reduction = "umap", group.by = "infercnv.group", 
        label = T, repel = T, pt.size = 0.5, label.size = 4.5) +
  labs(title = paste(sample1,"infercnv")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 13, face = "bold"))



integrated.2 <- integrated.2@meta.data %>%
  mutate(infercnv.group_sample = ifelse(paste(integrated.2@meta.data$Sample, 
                                               integrated.2@meta.data$infercnv.group,
                                               sep = "_")))


# Save the merged matrix of counts
#saveRDS(integrated.2, file = paste(prefix2, "merged_matrix_counts_cca_integration.rds", sep = "_")) 

#rm(anchors, anchors.2, integrated, integrated.2)
gc()


FeaturePlot(integrated.2, features = c("AGER","HMGB1", "S100A9", "S100A8", "S100A4", "S10012", "S100B"), 
            pt.size = 0.5, ncol = 3) &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

VlnPlot(integrated.2, features = c("AGER","HMGB1", "S100A9", "S100A8", "S100A4", "S100B"),
        ncol = 3, group.by = "cell_idents3", split.by = "Sample",
        cols = c("firebrick1", "hotpink", "purple", "cyan", "blue","green", "cadetblue4", "chartreuse4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

DotPlot(integrated.2, features = c("AGER","HMGB1", "S100A9", "S100A8", "S100A4", "S100B"),
        group.by = "cell_idents3",
        ) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

FeaturePlot(integrated.2, features = "S100A9", 
            pt.size = 0.5) +
  theme(axis.title = element_text(size = 9, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 9, face = "bold"))

FeaturePlot(integrated.2, features = "SPON1", 
            pt.size = 0.5) +
  theme(plot.title = element_text(size = 12))
  theme(axis.title = element_text(size = 9, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 9, face = "bold"))

v <- VlnPlot(integrated.2, features = "AGER",
       ncol = 1, group.by = "HPCA_labels_main") +
   theme(plot.title = element_text(size = 12)) +
   theme(axis.title = element_text(size = 9, face = "bold")) +
   theme(axis.line = element_line(size = 1.5)) +
   theme(axis.text = element_text(size=11, face = "bold")) +
   theme(legend.text = element_text(size = 10, face = "bold"))

v1 <- VlnPlot(integrated.2, features = "SPON1",
        ncol = 1, group.by = "HPCA_labels_main") +
   theme(plot.title = element_text(size = 12)) +
   theme(axis.title = element_text(size = 9, face = "bold")) +
   theme(axis.line = element_line(size = 1.5)) +
   theme(axis.text = element_text(size=11, face = "bold")) +
   theme(legend.text = element_text(size = 10, face = "bold"))

CombinePlots(plots = list(v, v1), legend = "right")

n <- FindMarkers(integrated.2, ident.1 = "Neurons",
                 ident.2 = c("Macrophages", "Monocytes",
                             "Neutrophils", "NK-cells", "T-cells",
                             "Tumor-Astrocytes", "Tumor-Oligodendrocytes"),
                 min.pct = 0.25, only.pos = T)
head(n, 40)


FeaturePlot(integrated.2, features = "NFASC")


head(integrated.2@meta.data, 5)


print("cell_type.txt file from annot")
head(integrated.2@meta.data, 5)
dim(integrated.2@meta.data)
celltype <- integrated.2@meta.data[,c(4,26)]
dim(celltype)
head(celltype,10)
setwd(dir = "D:/em_bild2/GBM/cca/")
#saveRDS(celltype, paste(prefix2, "cell_metadata.UMAPcluster.singleR_annot_cell_types.rds", sep="_"))
write.table(celltype, paste(prefix2, "cell_metadata.UMAPcluster.singleR_annot_cell_types.txt", sep="_"), 
            sep = "\t", quote = F, row.names = F, col.names = F)




annot_int <- integrated.2@meta.data

ggplot(as.data.frame(table(sort(annot_int$"Sample", decreasing = T))), aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill = "#1f1d66") +
  theme_classic(base_size = 18) + xlab("Sample") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 1.5)) +
  labs(title = "CELLS PER SAMPLE") +
  theme(plot.title = element_text(hjust = 0.5))


#correlation
gene_list <- c("AGER","HMGB1", "S100A9", "S100A8", "S100A4", "S100B",
               "IL1B", "TNF", "IL6", "GLO1", "IL10", "SPON1")
#correlation corrplot colors
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue")) 
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
whiteblack <- c("magenta", "blue")

#correlation integrated data
exp_mat <- FetchData(integrated.2, 
                       vars = c(gene_list),
                       slot = "data" )
#?cor.mtest
cor_exp_mat <- cor(exp_mat)
res1 <- cor.mtest(cor_exp_mat, conf.level = .95, method = "pearson")

corrplot.mixed(cor_exp_mat, lower.col = "black", order = "hclust",
               hclust.method = "ward.D",
               p.mat = res1$p, sig.level = 0.05, insig = "n",
               )

corrplot.mixed(cor_exp_mat, lower.col = "black", order = "hclust",
               hclust.method = "ward.D",
               p.mat = res1$p, sig.level = 0.05, insig = "label_sig",
)


FeatureScatter(integrated.2, feature1 = "S100A9", feature2 = "S100A8",
               group.by = "HPCA_labels_main", pt.size = 1,
               cols = c("firebrick1", "cadetblue1", "blue", "cadetblue4", "aquamarine","chartreuse4")) + 
  labs(title = paste(prefix2, "integation_AGER_vs_SPON1", sep = "_")) +
  labs(color = "Celltypes") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.text = element_text(size = 12, face = "bold"))


write.table(as.data.frame(table(sort(integrated.2$HPCA_labels_main))), 
            file = paste(prefix2, "celltype_numbers.txt", sep = '_'), 
            sep = "\t", quote = F)
####################################################################

#t-test

gene_list1 <- "AGER"

ager_test <- FetchData(integrated.2, vars = gene_list1,
                       slot = "data")
cell_names <- integrated.2@meta.data$HPCA_labels_main
ager_test[["Celltypes"]] <- cell_names
ager_test <- ager_test[, c("Celltypes","AGER")]
ager_test$Celltypes <- factor(ager_test$Celltypes, 
                              levels = c("Astrocyte", "Macrophage",
                                         "Monocyte", "NK_cell", 
                                         "Neutrophils", "T_cells"))
ager_test$Celltypes

baseformula <- ager_test$Celltypes
for (i in 2:ncol(ager_test)) {
  formula <- paste(colnames(ager_test) [i], 
                   baseformula, sep = "")
  
  p <- summary(aov(as.formula(formula), 
                   data = ager_test)) [[1]][["Pr(>F)"]][1]
  print(paste(formula, ":p=", p, sep = ""))
}

head(integrated.2@meta.data, 10)

for (i in 2:ncol(ager_test)) {
  print(colnames(ager_test)[i])
  print(summary(aov(ager_test[,i] ~ ager_test$Celltypes)))
}

gene_list1 <- "AGER"

ager_test <- FetchData(integrated.2, vars = gene_list1,
                       slot = "data")
cell_names <- integrated.2@meta.data$HPCA_labels_main
ager_test[["Celltypes"]] <- cell_names
ager_test <- ager_test[, c("Celltypes","AGER")]
ager_test <- within(ager_test, {
  Celltypes <- factor(Celltypes)
})

ager_test_mean <- aggregate(ager_test$AGER,
                            by = list(ager_test$Celltypes),
                            FUN = "mean")
colnames(ager_test_mean) <- c("Celltypes", "AGER")
head(ager_test_mean)
ager_aov <- with(ager_test_mean,
                 aov(AGER ~ Celltypes))
summary(ager_aov)

ag <- FindMarkers(integrated.2, test.use = "t",
                  ident.1 = "Monocytes", 
                  ident.2 = c("Malignant Astrocytes","Oligodendrocytes","Neurons",
                              "Macrophages","Neutrophils","T-cells","NK-cells"), 
                  features = "S100A9")

ag

write.table(ag, file = paste(prefix2, "S100A9_expression_monocytes_t.test.txt",
                             sep = "_"),
            sep = "\t", quote = F)
?FindMarkers
DimPlot(object = integrated.2, reduction = "umap", group.by = "RNA_snn_res.0.5", 
        label = T, pt.size = 0.5,
        cols = c("grey0", "grey0", "cadetblue1", "blue", "cadetblue4", "aquamarine","chartreuse4",
                 "darkslategray1", "darkorchid", "deeppink", "khaki1", "coral", "brown1", "mediumpurple")) +
  labs(title = paste(prefix2, "Seurat Clusters", sep = "_")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 14, face = "bold"))


myTheme <- theme_bw() +
  theme(text = element_text(size=16), axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "outside")


?DotPlot
DotPlot(integrated.2, features = c("AGER", "S100A9", "S100A8"))
head(integrated.2@meta.data, 5)



annot_int <- integrated.2@meta.data

setwd(dir = "D:/em_bild2/GBM/cca/figures/")

write.table(sort(table(annot_int$"Sample"), decreasing = T), 
            file = paste(prefix2, "cells_per_sample.txt", sep = "_"), 
            sep = "\t", quote = F, row.names = F)

ggplot(as.data.frame(table(annot_int$"Sample")), aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill = "#1f1d66") +
  theme_classic(base_size = 18) + xlab("Sample") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 1.5)) +
  labs(title = "CELLS PER SAMPLE") +
  theme(plot.title = element_text(hjust = 0.5)) +
  myTheme

pdf("001_makePreCCASeurat_Genes.pdf", width = 20, height = 6)
totalMedian <- median(integrated.2@meta.data$nFeature_RNA)
medianReads <- integrated.2@meta.data %>%
  group_by(Sample) %>%
  summarize(median = median(nFeature_RNA))
ggplot(medianReads, aes(x = Sample, y = median)) +
  ggtitle(paste(prefix2, "nFeature_RNA", sep = "_")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_bar(stat="identity", fill = "#1f1d66") + 
  geom_hline(aes(yintercept = totalMedian), lty = 2, col = "red") + 
  xlab("Sample name") + 
  ylab("Median distinct genes") + 
  myTheme
rm(medianReads, totalMedian)
dev.off()


pdf("001_makePreCCASeurat_Reads.pdf", width = 20, height = 6)
totalMedian <- median(integrated.2@meta.data$nCount_RNA)
medianReads <- integrated.2@meta.data %>%
  group_by(Sample) %>%
  summarize(median = median(nCount_RNA))
ggplot(medianReads, aes(x = Sample, y = median)) +
  ggtitle(paste(prefix2, "nCount_RNA", sep = "_")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_bar(stat="identity", fill = "#1f1d66") + 
  geom_hline(aes(yintercept = totalMedian), lty = 2, col = "red") + 
  xlab("Sample name") + 
  ylab("Median read coverage") + 
  myTheme

rm(medianReads, totalMedian)
dev.off()


pdf("001_makePreCCASeurat_Mitochondria.pdf", width = 20, height = 6)
ggplot(integrated.2@meta.data, aes(x = Sample, y = Percent.Mitochondria)) + 
  geom_boxplot() +
  ggtitle(paste(prefix2, "percent.mt", sep = "_")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Sample name") + 
  ylab("% Mitochondria") + 
  myTheme
dev.off()



head(integrated.2@meta.data, 5)


ma <- FindMarkers(integrated.2, ident.1 = "Malignant Astrocytes", test.use = "t",
                  ident.2 = c("Neurons", "Oligodendrocytes",
                              "Macrophages","Monocytes", "Neutrophils",
                              "NK-cells", "T-cells"),
                  min.pct = 0.5, only.pos = T)
head(ma, 30)



annot <- integrated.2@meta.data
dim(annot)
annot[["Sample2"]] <- integrated.2@meta.data$Sample
head(annot, 10)

annot <- within(annot,
                Sample2[Sample2 == "NP4396-Fluid1" & Sample == "NP4396-Fluid1"] <- "Fluid")
annot <- within(annot,
                Sample2[Sample2 == "NP4396-Fluid2" & Sample == "NP4396-Fluid2"] <- "Fluid")
annot <- within(annot,
                Sample2[Sample2 == "NP4396-TumorCenter" & Sample == "NP4396-TumorCenter"] <- "Tumor_Tissue")
annot <- within(annot,
                Sample2[Sample2 == "NP4396-TumorEdge" & Sample == "NP4396-TumorEdge"] <- "Tumor_Tissue")

annot <- within(annot,
                Sample2[Sample2 == "UR1113-Fluid1" & Sample == "UR1113-Fluid1"] <- "Fluid")
annot <- within(annot,
                Sample2[Sample2 == "UR1113-Fluid2" & Sample == "UR1113-Fluid2"] <- "Fluid")
annot <- within(annot,
                Sample2[Sample2 == "UR1113-TumorCenter" & Sample == "UR1113-TumorCenter"] <- "Tumor_Tissue")
annot <- within(annot,
                Sample2[Sample2 == "UR1113-TumorEdge" & Sample == "UR1113-TumorEdge"] <- "Tumor_Tissue")

annot <- within(annot,
                Sample2[Sample2 == "JS4587_fluid" & Sample == "JS4587_fluid"] <- "Fluid")
annot <- within(annot,
                Sample2[Sample2 == "JS4587_tumor_center" & Sample == "JS4587_tumor_center"] <- "Tumor_Tissue")
annot <- within(annot,
                Sample2[Sample2 == "JS4587_tumor_edge" & Sample == "JS4587_tumor_edge"] <- "Tumor_Tissue")

annot <- within(annot,
                Sample2[Sample2 == "RY8558_fluid" & Sample == "RY8558_fluid"] <- "Fluid")
annot <- within(annot,
                Sample2[Sample2 == "RY8558_tumor_center" & Sample == "RY8558_tumor_center"] <- "Tumor_Tissue")
annot <- within(annot,
                Sample2[Sample2 == "RY8558_tumor_edge" & Sample == "RY8558_tumor_edge"] <- "Tumor_Tissue")

annot <- within(annot,
                Sample2[Sample2 == "YY2224_fluid" & Sample == "YY2224_fluid"] <- "Fluid")
annot <- within(annot,
                Sample2[Sample2 == "YY2224_tumor_center" & Sample == "YY2224_tumor_center"] <- "Tumor_Tissue")
annot <- within(annot,
                Sample2[Sample2 == "YY2224_tumor_edge" & Sample == "YY2224_tumor_edge"] <- "Tumor_Tissue")

table(sort(annot$Sample2, decreasing = T))

head(annot, 10)
dim(annot)
head(integrated.2@meta.data, 10)
dim(integrated.2@meta.data)

integrated.2[["Sample2"]] <- annot$Sample2

VlnPlot(integrated.2, features = c("AGER","HMGB1", "S100A9", "S100A8", "S100A4", "S100B"),
        ncol = 3, group.by = "cell_idents3", split.by = "Sample2",
        cols = c("darkslategray1","green4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold")) 

integrated.2_Monocytes <- subset(integrated.2, 
                                 subset = cell_idents3 != "Malignant Astrocytes" )
integrated.2_Monocytes <- subset(integrated.2_Monocytes, 
                                 subset = cell_idents3 != "Oligodendrocytes" )
integrated.2_Monocytes <- subset(integrated.2_Monocytes, 
                                 subset = cell_idents3 != "Neurons" )
integrated.2_Monocytes <- subset(integrated.2_Monocytes, 
                                 subset = cell_idents3 != "Macrophages" )
integrated.2_Monocytes <- subset(integrated.2_Monocytes, 
                                 subset = cell_idents3 != "Neutrophils" )
integrated.2_Monocytes <- subset(integrated.2_Monocytes, 
                                 subset = cell_idents3 != "T-cells" )
integrated.2_Monocytes <- subset(integrated.2_Monocytes, 
                                 subset = cell_idents3 != "NK-cells" )

table(sort(integrated.2_Monocytes$cell_idents3))
dim(integrated.2_Monocytes)
integrated.2_Monocytes$celltype.condition <- paste(Idents(integrated.2_Monocytes), 
                                                   integrated.2_Monocytes$Sample2, sep="_")
integrated.2_Monocytes$celltype <- Idents(integrated.2_Monocytes)
Idents(integrated.2_Monocytes) <- "celltype.condition"
Idents(integrated.2_Monocytes)

#for (i in 0:2){ #or however many clusters you have
#  try({
#    ident1 <- paste0(i,"_condition1")
#    ident2 <- paste0(i,"condition2")
#    condition.diffgenes <- FindMarkers(integrated.2_Monocytes, ident.1 = ident1, ident.2=ident2, min.pct=0.25, logfc.threshold=0.25)
#  })
#}

ag <- FindMarkers(integrated.2_Monocytes, test.use = "t",
                  ident.1 = "Monocytes_Fluid", 
                  ident.2 = c("Monocytes_Tumor_Tissue"), 
                  features = "S100A9")

ag


VlnPlot(integrated.2_Monocytes, features = c("AGER","HMGB1", "S100A9", "S100A8", "S100A4", "S100B"),
        ncol = 3,
        cols = c("darkslategray1","green4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

VlnPlot(integrated.2_Monocytes, features = "S100A9",
        cols = c("darkslategray1","green4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))


setwd(dir = "D:/em_bild2/GBM/cca/seurat_objects3/")
#saveRDS(integrated.2_Monocytes, paste(prefix2, "cca_seurat_obj_Monocytes_only.rds", sep = "_"))
integrated.2_Monocytes <- readRDS(paste(prefix2, "cca_seurat_obj_Monocytes_only.rds", sep = "_"))

table(sort(integrated.2_Monocytes$Sample2, decreasing = T))
ag <- FindMarkers(integrated.2_fluid, test.use = "t",
                  ident.1 = "Monocytes", 
                  ident.2 = c("Malignant Astrocytes","Oligodendrocytes",
                              "Macrophages","Neutrophils","T-cells","NK-cells"), 
                  features = "S100A9")

ag


head(integrated.2_Monocytes@meta.data, 10)

annot <- integrated.2_Monocytes@meta.data
head(annot, 5)
annot[["Sample3"]] <- annot$Sample
table(sort(annot$Sample, decreasing = T))
annot <- within(annot,
                Sample3[Sample3 == "JS4587_fluid" & Sample == "JS4587_fluid"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "RY8558_fluid" & Sample == "RY8558_fluid"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "YY2224_fluid" & Sample == "YY2224_fluid"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "NP4396-Fluid1" & Sample == "NP4396-Fluid1"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "UR1113-Fluid1" & Sample == "UR1113-Fluid1"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "NP4396-Fluid2" & Sample == "NP4396-Fluid2"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "UR1113-Fluid2" & Sample == "UR1113-Fluid2"] <- "Fluid")

annot <- within(annot,
                Sample3[Sample3 == "JS4587_tumor_center" & Sample == "JS4587_tumor_center"] <- "Tumor_Center")
annot <- within(annot,
                Sample3[Sample3 == "RY8558_tumor_center" & Sample == "RY8558_tumor_center"] <- "Tumor_Center")
annot <- within(annot,
                Sample3[Sample3 == "YY2224_tumor_center" & Sample == "YY2224_tumor_center"] <- "Tumor_Center")
annot <- within(annot,
                Sample3[Sample3 == "NP4396-TumorCenter" & Sample == "NP4396-TumorCenter"] <- "Tumor_Center")
annot <- within(annot,
                Sample3[Sample3 == "UR1113-TumorCenter" & Sample == "UR1113-TumorCenter"] <- "Tumor_Center")


annot <- within(annot,
                Sample3[Sample3 == "JS4587_tumor_edge" & Sample == "JS4587_tumor_edge"] <- "Tumor_Edge")
annot <- within(annot,
                Sample3[Sample3 == "RY8558_tumor_edge" & Sample == "RY8558_tumor_edge"] <- "Tumor_Edge")
annot <- within(annot,
                Sample3[Sample3 == "YY2224_tumor_edge" & Sample == "YY2224_tumor_edge"] <- "Tumor_Edge")
annot <- within(annot,
                Sample3[Sample3 == "NP4396-TumorEdge" & Sample == "NP4396-TumorEdge"] <- "Tumor_Edge")
annot <- within(annot,
                Sample3[Sample3 == "UR1113-TumorEdge" & Sample == "UR1113-TumorEdge"] <- "Tumor_Edge")

table(sort(annot$Sample3, decreasing = T))

integrated.2_Monocytes[["Sample3"]] <- annot$Sample3
table(sort(integrated.2_Monocytes$Sample3, decreasing = T))


table(sort(integrated.2$Sample))
rm(integrated.2_tumor)
integrated.2_fluid <- subset(integrated.2, 
                             subset = Sample != "JS4587_tumor_center")
integrated.2_fluid <- subset(integrated.2_fluid, 
                             subset = Sample != "JS4587_tumor_edge")                             
integrated.2_fluid <- subset(integrated.2_fluid, 
                             subset = Sample != "RY8558_tumor_center")
integrated.2_fluid <- subset(integrated.2_fluid, 
                             subset = Sample != "RY8558_tumor_edge") 
integrated.2_fluid <- subset(integrated.2_fluid, 
                             subset = Sample != "YY2224_tumor_center")
integrated.2_fluid <- subset(integrated.2_fluid, 
                             subset = Sample != "YY2224_tumor_edge")
integrated.2_fluid <- subset(integrated.2_fluid, 
                             subset = Sample != "NP4396-TumorCenter")
integrated.2_fluid <- subset(integrated.2_fluid, 
                             subset = Sample != "NP4396-TumorEdge")
integrated.2_fluid <- subset(integrated.2_fluid, 
                             subset = Sample != "UR1113-TumorCenter")
integrated.2_fluid <- subset(integrated.2_fluid, 
                             subset = Sample != "UR1113-TumorEdge")

integrated.2_fluid
table(sort(integrated.2_fluid$Sample))
setwd(dir = "/Volumes/T7/em_bild2/GBM/cca/seurat_objects3/")
#saveRDS(integrated.2_fluid, paste(prefix2, "cca_seurat_obj_fluid_only.rds", sep = "_"))
integrated.2_fluid <- readRDS(file = paste(prefix2, "cca_seurat_obj_fluid_only.rds", sep = "_"))
integrated.2_fluid


VlnPlot(integrated.2_fluid, features = c("AGER","HMGB1", "S100A9", "S100A8", "S100A4", "S100B"),
        ncol = 3, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "cyan", "blue","green", "cadetblue4", "chartreuse4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

ag <- FindMarkers(integrated.2_fluid, test.use = "t",
                  ident.1 = "Monocytes", 
                  ident.2 = c("Malignant Astrocytes","Oligodendrocytes",
                              "Macrophages","Neutrophils","T-cells","NK-cells"), 
                  features = "S100A9")

ag


VlnPlot(integrated.2, features = c("LGALS3"),
        ncol = 1, group.by = "cell_idents3", split.by = "Sample3",
        cols = c("darkslategray1","green4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))



table(sort(integrated.2$Sample))

integrated.2_tumor <- subset(integrated.2, 
                             subset = Sample != "JS4587_fluid")
integrated.2_tumor <- subset(integrated.2_tumor, 
                             subset = Sample != "RY8558_fluid")
integrated.2_tumor <- subset(integrated.2_tumor, 
                             subset = Sample != "YY2224_fluid")
integrated.2_tumor <- subset(integrated.2_tumor, 
                             subset = Sample != "NP4396-Fluid1")
integrated.2_tumor <- subset(integrated.2_tumor, 
                             subset = Sample != "NP4396-Fluid2")
integrated.2_tumor <- subset(integrated.2_tumor, 
                             subset = Sample != "UR1113-Fluid1")
integrated.2_tumor <- subset(integrated.2_tumor, 
                             subset = Sample != "UR1113-Fluid2")
integrated.2_tumor
table(sort(integrated.2_tumor$Sample))

setwd(dir = "/Volumes/T7//em_bild2/GBM/cca/seurat_objects3/")
#saveRDS(integrated.2_tumor, paste(prefix2, "cca_seurat_obj_tumor_only.rds", sep = "_"))
integrated.2_tumor <- readRDS(file = paste(prefix2, "cca_seurat_obj_tumor_only.rds", sep = "_"))
integrated.2_tumor


VlnPlot(integrated.2_tumor, features = c("AGER","HMGB1", "S100A9", "S100A8", "S100A4", "S100B"),
        ncol = 3, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "purple", "cyan", "blue","green", "cadetblue4", "chartreuse4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

###individual genes
VlnPlot(integrated.2_tumor, features = "S100A9",
        ncol = 1, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "purple", "cyan", "blue","green", "cadetblue4", "chartreuse4")) +
  labs(x="") +
  theme(axis.title = element_text(size = 9, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 9, face = "bold")) +
  NoLegend()
VlnPlot(integrated.2_tumor, features = "HMGB1",
        ncol = 1, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "purple", "cyan", "blue","green", "cadetblue4", "chartreuse4")) +
  labs(x="") +
  theme(axis.title = element_text(size = 9, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 9, face = "bold")) +
  NoLegend()
VlnPlot(integrated.2_tumor, features = "S100B",
        ncol = 1, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "purple", "cyan", "blue","green", "cadetblue4", "chartreuse4")) +
  labs(x="") +
  theme(axis.title = element_text(size = 9, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 9, face = "bold")) +
  NoLegend()

ag <- FindMarkers(integrated.2_tumor, test.use = "t",
                  ident.1 = "Monocytes", 
                  ident.2 = c("Malignant Astrocytes","Oligodendrocytes","Neurons",
                              "Macrophages","Neutrophils","T-cells","NK-cells"), 
                  features = "S100A9")

ag


FeaturePlot(integrated.2, features = "LGALS3")

annot <- integrated.2@meta.data
dim(annot)
annot[["Sample3"]] <- integrated.2@meta.data$Sample
head(annot, 10)

annot <- within(annot,
                Sample3[Sample3 == "NP4396-Fluid1" & Sample == "NP4396-Fluid1"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "NP4396-Fluid2" & Sample == "NP4396-Fluid2"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "NP4396-TumorCenter" & Sample == "NP4396-TumorCenter"] <- "Tumor_Center")
annot <- within(annot,
                Sample3[Sample3 == "NP4396-TumorEdge" & Sample == "NP4396-TumorEdge"] <- "Tumor_Edge")

annot <- within(annot,
                Sample3[Sample3 == "UR1113-Fluid1" & Sample == "UR1113-Fluid1"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "UR1113-Fluid2" & Sample == "UR1113-Fluid2"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "UR1113-TumorCenter" & Sample == "UR1113-TumorCenter"] <- "Tumor_Center")
annot <- within(annot,
                Sample3[Sample3 == "UR1113-TumorEdge" & Sample == "UR1113-TumorEdge"] <- "Tumor_Edge")

annot <- within(annot,
                Sample3[Sample3 == "JS4587_fluid" & Sample == "JS4587_fluid"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "JS4587_tumor_center" & Sample == "JS4587_tumor_center"] <- "Tumor_Center")
annot <- within(annot,
                Sample3[Sample3 == "JS4587_tumor_edge" & Sample == "JS4587_tumor_edge"] <- "Tumor_Edge")

annot <- within(annot,
                Sample3[Sample3 == "RY8558_fluid" & Sample == "RY8558_fluid"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "RY8558_tumor_center" & Sample == "RY8558_tumor_center"] <- "Tumor_Center")
annot <- within(annot,
                Sample3[Sample3 == "RY8558_tumor_edge" & Sample == "RY8558_tumor_edge"] <- "Tumor_Edge")

annot <- within(annot,
                Sample3[Sample3 == "YY2224_fluid" & Sample == "YY2224_fluid"] <- "Fluid")
annot <- within(annot,
                Sample3[Sample3 == "YY2224_tumor_center" & Sample == "YY2224_tumor_center"] <- "Tumor_Center")
annot <- within(annot,
                Sample3[Sample3 == "YY2224_tumor_edge" & Sample == "YY2224_tumor_edge"] <- "Tumor_Edge")

table(sort(annot$Sample3, decreasing = T))

head(annot, 10)
dim(annot)
head(integrated.2@meta.data, 10)
dim(integrated.2@meta.data)

integrated.2[["Sample3"]] <- annot$Sample3


FeaturePlot(integrated.2, features = "LGALS3")


VlnPlot(integrated.2, features = c("LGALS3"),
        ncol = 1, group.by = "cell_idents3", split.by = "Sample3",
        cols = c("darkslategray1","darkorchid", "deeppink")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold")) 

VlnPlot(integrated.2_fluid, features = c("LGALS3"),
        ncol = 1, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "cyan", "blue","green", "cadetblue4", "chartreuse4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

FeaturePlot(integrated.2, features = "LGALS3", 
            pt.size = 0.5) +
  theme(axis.title = element_text(size = 9, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 9, face = "bold"))



###individual genes in fluid only
VlnPlot(integrated.2_fluid, features = "S100A9",
        ncol = 1, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "cyan", "blue","green", "cadetblue4", "chartreuse4")) +
  labs(x="") +
  theme(axis.title = element_text(size = 9, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 9, face = "bold")) +
  NoLegend()
VlnPlot(integrated.2_fluid, features = "HMGB1",
        ncol = 1, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "cyan", "blue","green", "cadetblue4", "chartreuse4")) +
  labs(x="") +
  theme(axis.title = element_text(size = 9, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 9, face = "bold")) +
  NoLegend()
VlnPlot(integrated.2_fluid, features = "S100B",
        ncol = 1, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "cyan", "blue","green", "cadetblue4", "chartreuse4")) +
  labs(x="") +
  theme(axis.title = element_text(size = 9, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 9, face = "bold")) +
  NoLegend()

table(sort(integrated.2_tumor$Sample, decreasing = T))
integrated.2_edge <- subset(integrated.2_tumor, 
                             subset = Sample != "JS4587_tumor_center")
integrated.2_edge <- subset(integrated.2_edge, 
                             subset = Sample != "RY8558_tumor_center")
integrated.2_edge <- subset(integrated.2_edge, 
                             subset = Sample != "YY2224_tumor_center")
integrated.2_edge <- subset(integrated.2_edge, 
                             subset = Sample != "NP4396-TumorCenter")
integrated.2_edge <- subset(integrated.2_edge, 
                             subset = Sample != "UR1113-TumorCenter")
table(sort(integrated.2_edge$Sample, decreasing = T))
setwd(dir = "D:/em_bild2/GBM/cca/seurat_objects3/")
saveRDS(integrated.2_edge, file = paste(prefix2, "cca_seurat_obj_tumor.edge_only.rds", sep = "_"))

VlnPlot(integrated.2_edge, features = c("LGALS3"),
        ncol = 1, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "purple", "cyan", "blue","green", "cadetblue4", "chartreuse4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))



integrated.2_center <- subset(integrated.2_tumor, 
                              subset = Sample != "JS4587_tumor_edge")  
integrated.2_center <- subset(integrated.2_center, 
                             subset = Sample != "RY8558_tumor_edge")
integrated.2_center <- subset(integrated.2_center, 
                             subset = Sample != "YY2224_tumor_edge")
integrated.2_center <- subset(integrated.2_center, 
                             subset = Sample != "NP4396-TumorEdge")
integrated.2_center <- subset(integrated.2_center, 
                             subset = Sample != "UR1113-TumorEdge")
table(sort(integrated.2_center$Sample, decreasing = T))
setwd(dir = "D:/em_bild2/GBM/cca/seurat_objects3/")
saveRDS(integrated.2_center, file = paste(prefix2, "cca_seurat_obj_tumor.center_only.rds", sep = "_"))
VlnPlot(integrated.2_center, features = c("LGALS3"),
        ncol = 1, group.by = "cell_idents3",
        cols = c("firebrick1", "hotpink", "purple", "cyan", "blue","green", "cadetblue4", "chartreuse4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))




head(integrated.2@meta.data, 10)

table(sort(annot$Sample3, decreasing = T))
table(sort(integrated.2$Sample3, decreasing = T))


VlnPlot(integrated.2, features = c("LGALS3"),
        ncol = 1, group.by = "cell_idents3", split.by = "Sample3",
        cols = c("firebrick1", "hotpink", "purple", "cyan", "blue","green", "cadetblue4", "chartreuse4")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))


ag <- FindMarkers(integrated.2, test.use = "t", group.by = "Sample3",
                  ident.1 = "Monocytes", 
                  ident.2 = c("Malignant Astrocytes","Oligodendrocytes","Neurons",
                              "Macrophages","Neutrophils","T-cells","NK-cells"), 
                  features = "LGALS3")

ag


table(sort(integrated.2_Monocytes$Sam))

VlnPlot(integrated.2_Monocytes, features = c("LGALS3"),
        ncol = 1, group.by = "cell_idents3", split.by = "Sample3",
        cols = c("darkslategray1","darkorchid", "deeppink")) &
  labs(x="") &
  theme(axis.title = element_text(size = 9, face = "bold")) &
  theme(axis.line = element_line(size = 1.5)) &
  theme(axis.text = element_text(size=10, face = "bold")) &
  theme(legend.text = element_text(size = 9, face = "bold"))

ag <- FindMarkers(integrated.2_Monocytes, test.use = "t", group.by = "Sample3",
                  ident.1 = "Fluid", 
                  ident.2 = c("Tumor_Center"), 
                  features = "LGALS3")

ag

ag <- FindMarkers(integrated.2_Monocytes, test.use = "t", group.by = "Sample3",
                  ident.1 = "Fluid", 
                  ident.2 = c("Tumor_Edge"), 
                  features = "LGALS3")

ag

annots <- integrated.2@meta.data
dim(annots)
np4396_metadata <- annots[grepl("NP4396", rownames(annots)),]

np4396_annotations <- matrix(nrow = 3131, ncol = 2)
np4396_annotations <- np4396_metadata[,c(4,27)]
np4396_annotations <- np4396_annotations %>% 
  mutate(cell_idents3 = ifelse(cell_idents3 == "Malignant Astrocytes", 
                               "Malignant Astrocytes", "normal"))

write.table(np4396_annotations, file = paste(prefix2, 
                                             "np4396_annotations_infercnv.txt",
                                             sep = "_"),
            sep = "\t", row.names = F, col.names = F, quote = F)
#
ur1113_metadata <- annots[grepl("UR1113", rownames(annots)),]
ur1113_annotations <- ur1113_metadata[,c(4,27)]
ur1113_annotations <- ur1113_annotations %>% 
  mutate(cell_idents3 = ifelse(cell_idents3 == "Malignant Astrocytes", 
                               "Malignant Astrocytes", "normal"))
write.table(ur1113_annotations, file = paste(prefix2, 
                                             "ur1113_annotations_infercnv.txt",
                                             sep = "_"),
            sep = "\t", row.names = F, col.names = F, quote = F)
#
js4587_metadata <- annots[grepl("JS4587", rownames(annots)),]
js4587_annotations <- js4587_metadata[,c(4,27)]
js4587_annotations <- js4587_annotations %>% 
  mutate(cell_idents3 = ifelse(cell_idents3 == "Malignant Astrocytes", 
                               "Malignant Astrocytes", "normal"))
write.table(js4587_annotations, file = paste(prefix2, 
                                             "js4587_annotations_infercnv.txt",
                                             sep = "_"),
            sep = "\t", row.names = F, col.names = F, quote = F)
ry8558_metadata <- annots[grepl("RY8558", rownames(annots)),]
ry8558_annotations <- ry8558_metadata[,c(4,27)]
ry8558_annotations <- ry8558_annotations %>% 
  mutate(cell_idents3 = ifelse(cell_idents3 == "Malignant Astrocytes", 
                               "Malignant Astrocytes", "normal"))
write.table(ry8558_annotations, file = paste(prefix2, 
                                             "ry8558_annotations_infercnv.txt",
                                             sep = "_"),
            sep = "\t", row.names = F, col.names = F, quote = F)
#
yy2224_metadata <- annots[grepl("YY2224", rownames(annots)),]
yy2224_annotations <- yy2224_metadata[,c(4,27)]
yy2224_annotations <- yy2224_annotations %>% 
  mutate(cell_idents3 = ifelse(cell_idents3 == "Malignant Astrocytes", 
                               "Malignant Astrocytes", "normal"))
write.table(yy2224_annotations, file = paste(prefix2, 
                                             "yy2224_annotations_infercnv.txt",
                                             sep = "_"),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(as.matrix(GetAssayData(object = integrated.2, slot = "counts")),
            file = paste(prefix2, "raw_spars_matrix_counts.txt", sep = "_"),
            sep = "\t", row.names = T, col.names = T, quote = F)
