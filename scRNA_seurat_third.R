print("Load R libraries")
library(dplyr)
library(data.table)
library(Seurat)
library(patchwork)
library(ggplot2)
library(utils)
library(RColorBrewer)
library(Matrix)

#detach("package:dplyr", unload = T)

if (FALSE) {stop("The value is TRUE, so the script must end here")
} else { #continue the script
  print("Script did NOT end!")
}

sample <- "YY2224"
sample0 <- "COH084"
sample1 <- paste(sample0, sample, sep = "_")
platform <- "10x"
prefix0 <- paste(sample0, platform, sep = "_")
prefix1 <- paste(sample1, platform, sep = "_")

print("load raw seurat obj from first run")
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat"))
seu1 <- readRDS(file = paste(prefix1,
                             "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run2_infercnv.rds", 
                             sep = "_"))

cl1 <- FindMarkers(seu1, ident.1 = "clone2", test.use = "t",
                   ident.2 = c("Macrophages", "Monocytes","T-cells"),
                   min.pct = 0.25, only.pos = T)
head(cl1, 40)


FindMarkers(seu1, ident.1 = "clone1", test.use = "t",
            ident.2 = c("Macrophages-1", "Macrophages-2","Monocytes","T-cells", "NK-cells"),
            features = "")


FeaturePlot(seu1, features = c("CD14", "AIF1", "FCER1G", "FCGR3A", "TYROBP", "CSF1R"))
FeaturePlot(seu1, features = c("MBP", "TF", "PLP1"))


genes <- matrix(nrow = 9, ncol = 2)
colnames(genes) <- c("Celltype", "Markers")
genes[1,1] <- "Tumor Astrocytes"
genes[1,2] <- c("GFAP,ALDH1L1")
genes[2,1] <- "Tumor Oligodendrocytes"
genes[2,2] <- c("PLP1,TF")
genes[3,1] <- "Neurons"
genes[3,2] <- c("SYT1")
genes[4,1] <- "Monocytes"
genes[4,2] <- c("S100A9, LYZ")
genes[3,1] <- "Macrophages"
genes[3,2] <- c("")
genes[4,1] <- "T-cells"
genes[4,2] <- c("")
genes[3,1] <- "NK-cells"
genes[3,2] <- c("")
genes[4,1] <- "Neutrophils"
genes[4,2] <- c("")

genes


genes <- matrix(nrow = 10, ncol = 2)
colnames(genes) <- c("Celltype", "Markers")
genes[1,1] <- "Tumor Astrocytes"
genes[1,2] <- c("GFAP,EGFR")
genes[2,1] <- "Tumor Oligodendrocytes"
genes[2,2] <- c("PLP1,TF")
genes[3,1] <- "Neurons"
genes[3,2] <- c("SYT1")
genes[4,1] <- "Monocytes"
genes[4,2] <- c("S100A9, LYZ")
genes[3,1] <- "Macrophages"
genes[3,2] <- c("")
genes[4,1] <- "T-cells"
genes[4,2] <- c("")
genes[3,1] <- "NK-cells"
genes[3,2] <- c("")
genes[4,1] <- "Neutrophils"
genes[4,2] <- c("")


genes <- matrix(nrow = 10, ncol = 2)
colnames(genes) <- c("Celltype", "Markers")
genes[1,1] <- "Tumor Astrocytes"
genes[1,2] <- c("GFAP")
genes[2,1] <- "Tumor Oligodendrocytes"
genes[2,2] <- c("PLP1,TF")
genes[3,1] <- "Neurons"
genes[3,2] <- c("SYT1")
genes[4,1] <- "Monocytes"
genes[4,2] <- c("S100A9, LYZ")
genes[3,1] <- "Macrophages"
genes[3,2] <- c("")
genes[4,1] <- "T-cells"
genes[4,2] <- c("")
genes[3,1] <- "NK-cells"
genes[3,2] <- c("")
genes[4,1] <- "Neutrophils"
genes[4,2] <- c("")

head(Idents(seu1),5)
seu2 <- seu1
seu2@meta.data$cell_idents <- as.character(seu2@meta.data$cell_idents)
seu2@meta.data <- within(seu2@meta.data, 
                        cell_idents[cell_idents == "Monocytes" & RNA_snn_res.0.5 == "8"] <- "Macrophages")
seu2@meta.data$cell_idents <- as.factor(seu2@meta.data$cell_idents)

Idents(seu2) <- seu2@meta.data$cell_idents
seu2 <- RenameIdents(object = seu2,
                     'Tumor-Astrocytes' = "Malignant-Astrocytes",
                     'Tumor-Oligodendrocytes' = "Oligodendrocytes",
                     'Macrophages' = "Macrophages",
                     'Macrophages-2' = "Macrophages",
                     'Monocytes' = "Monocytes",
                     'T-cells' = "T-cells",
                     'NK-cells' = "NK-cells",
                     'clone2' = "Neurons",
                     'Neutrophils' = "Neutrophils")

seu2[["cell_idents3"]] <- Idents(object = seu2)
head(Idents(seu2), 5)

DimPlot(seu2, reduction = "umap", group.by = "cell_idents3", label = F, 
        repel = T, pt.size = 0.5, label.size = 3.5) +
  labs(title = paste(prefix1, "Cell Idents", sep = " ")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 13, face = "bold"))

DimPlot(seu2, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T, 
        repel = T, pt.size = 0.5, label.size = 3.5) +
  labs(title = paste(prefix1, "Clusters", sep = " ")) +
  theme(axis.title = element_text(size = 10, face = "bold")) +
  theme(axis.line = element_line(size = 1.5)) +
  theme(axis.text = element_text(size=10, face = "bold")) +
  theme(legend.text = element_text(size = 13, face = "bold"))

seu1 <- seu2
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat"))
table(sort(seu1@meta.data$cell_idents3, decreasing = T))
print("cell_type.txt file from annot")
head(seu1@meta.data, 5)
dim(seu1@meta.data)
celltype <- seu1@meta.data[,c(4,26)]
dim(celltype)
head(celltype,10)
saveRDS(celltype, paste(prefix1, "cell_metadata.UMAPcluster.singleR_annot_cell_types3.rds", sep="_"))
write.table(celltype, paste(prefix1, "cell_metadata.UMAPcluster.singleR_annot_cell_types3.txt", sep="_"), 
            sep = "\t", quote = F, row.names = F, col.names = F)
#saveRDS(seu1, file = paste(prefix1, "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run3_infercnv.rds", sep="_"))


seu1 <- readRDS(file = paste(prefix1,
                             "Seurat_2kgenes_vst_cc.final.seu.obj_changed_idents_sec_run3_infercnv.rds", 
                             sep = "_"))


genes <- matrix(nrow = 8, ncol = 2)
colnames(genes) <- c("Celltype", "Markers")
genes[1,1] <- "Malignant Astrocytes"
genes[1,2] <- c("GFAP")
genes[2,1] <- "Tumor Oligodendrocytes"
genes[2,2] <- c("PLP1,TF")
genes[3,1] <- "Neurons"
genes[3,2] <- c("SYT1")
genes[4,1] <- "Monocytes"
genes[4,2] <- c("S100A9, LYZ")
genes[3,1] <- "Macrophages"
genes[3,2] <- c("")
genes[4,1] <- "T-cells"
genes[4,2] <- c("")
genes[3,1] <- "NK-cells"
genes[3,2] <- c("")
genes[4,1] <- "Neutrophils"
genes[4,2] <- c("")

genes

ma <- FindMarkers(integrated.2, ident.1 = "Malignant Astrocytes",
                  ident.2 = c(""),
                  min.pct = 0.5, only.pos = T)
head(ma, 10)

print("END OF SCRIPT")

