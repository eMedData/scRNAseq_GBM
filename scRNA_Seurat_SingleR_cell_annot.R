print("loading libraries")
library(SingleR)
library(scRNAseq)
library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(utils)
library(viridis)
library(celldex)

detach("package:", unload=TRUE)

#switch out prefix1 samples in .sh file before running


sample <- "YY2224"
sample0 <- "COH084"
sample1 <- paste(sample0, sample, sep = "_")
platform <- "10x"
prefix0 <- paste(sample0, platform, sep = "_")
prefix1 <- paste(sample1, platform, sep = "_")

print("loading counts")
#use this function to skip chunk of code
.f = function() {
counts.df <-  fread(paste(sample, platform, "count01.noaggr.counts.txt", sep="_"), sep="\t")
setwd(dir = "/net/isi-dcnl/ifs/user_data/abild/eric/COH086/")
counts.df <- fread("COH086_10x_read_counts.txt", sep = "\t")
counts.df <- subset(counts.df, counts.df$'Gene Symbol' != "")
counts.df <- unique(counts.df, by='Gene Symbol')
counts.df <- as.data.frame(counts.df)
rownames(counts.df) <- counts.df$`Gene Symbol`
counts <- counts.df[,4:ncol(counts.df)]
counts[1:5,1:5]
}

setwd(dir = "D:/em_bild2/GBM/")
counts <- readRDS(file = paste(prefix0, "gene.symbol_read_counts.rds", sep="_"))
dim(counts)
counts[1:5, 1:5]
print("done")

#Annotation
print("loading cellmetadata")
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/seurat"))
annot <- readRDS(paste(prefix1, "cell_metadata.UMAPcluster_first_run.rds", sep="_"))
dim(annot)
rownames(annot) <- annot$Cell
annot[1:5,]
print("done")

print("matching counts to cellmetadata")
dim(counts)
counts <- counts[match(annot$"Cell", colnames(counts))]
dim(counts)
counts[1:5,1:5]

print("matching cellmetadata to counts")
dim(annot)
annot <- annot[match(colnames(counts), annot$"Cell"),]
dim(annot)
annot[1:5,]

print("loading singleR refs")
#ref.hpca is reference data set (Mabbott et al. 2013)
ref.hpca <- celldex::HumanPrimaryCellAtlasData()
#cell types
colnames(ref.hpca[,1:10])
#genes
rownames(ref.hpca[1:10,])
print("done")

#ref.encode (BlueprintEncodeData)
ref.encode <- celldex::BlueprintEncodeData(rm.NA = "rows")
#cell types
colnames(ref.encode[,1:10])
#genes
rownames(ref.encode[1:10,])


print("counts as matrix")
typeof(counts)
counts <- as.matrix(counts)
typeof(counts)
print("done")


if (FALSE) {stop("The value is TRUE, so the script must end here")
  } else { #continue the script
print("Script did NOT end!")
  }

print("singleR using hpca main types start")
table(sort(ref.hpca$label.main, decreasing = T))
singler <- SingleR(test = counts, ref = ref.hpca, labels = ref.hpca$label.main, 
                   clusters = NULL, genes = "sd", sd.thresh = 1, 
                   quantile = 0.8, fine.tune = T, tune.thresh = 0.05)


print("saving hpca.main singleR obj")
setwd(dir = paste0("D:/em_bild2/GBM/", prefix1, "/singleR"))
saveRDS(singler, file = paste(prefix1, "singleR.object_hpca_labels_main.rds", sep="_"))
sort(table(singler$labels), decreasing=T)
singler[1:10, 1:5]
write.table(sort(table(singler$labels), decreasing = T), 
            file = paste(prefix1, "singleR_ref.hpca_labels_main_table.txt", sep="_"), 
            sep="\t", quote = F, row.names = F)



print("saving heatmaps")
plotScoreHeatmap(singler, clusters = annot$Sample, fontsize = 12)
#plotScoreHeatmap(singler, clusters = annot$RNA_snn_res.0.5, fontsize = 12)
#plotDeltaDistribution(singler, ncol = 4)
print("done")

print("inputting cell types into metadata")
y <- match(colnames(counts), rownames(singler))
annot[["HPCA_labels_main"]] <- singler$labels[y]
dim(annot)
annot[1:5,]
print("singleR with hpca main types done")


print("start encode ref types")
table(sort(ref.encode$label.main, decreasing = T))
singler <- SingleR(test = counts, ref = ref.encode, labels = ref.encode$label.main, 
                   genes = "sd", sd.thresh = 1, quantile = 0.8, 
                   fine.tune = T, tune.thresh = 0.05)

print("saving encode types singleR obj")
saveRDS(singler, file = paste(prefix1, "singleR.object_encode_labels_main.rds", sep="_"))
write.table(sort(table(singler$labels), decreasing = T), 
            file = paste(prefix1, "singleR_encode_labels_main_table.txt", sep="_"), sep="\t", quote = F, row.names = F)

print("saving heatmaps")
plotScoreHeatmap(singler, clusters = annot$Sample, fontsize = 12)

print("inputting cell types into metadata")
singler[1:5, 1:5]
sort(table(singler$labels), decreasing=T)
y <- match(colnames(counts), rownames(singler))
annot[["ENCODE_labels_main"]] <- singler$labels[y]
dim(annot)
annot[1:5,]

print("saving _cell_annot.txt file")
write.table(annot, paste(prefix1, "cell_metadata.UMAPcluster.annot_singleR.txt", sep="_"), 
            sep = "\t", quote = F, row.names= F)
print("cell_type.txt file from annot")
celltype <- annot[,c(1,6,21,22)]
celltype
saveRDS(celltype, paste(prefix1, "cell_metadata.UMAPcluster.singleR_annot_cell_types.rds", sep="_"))
write.table(celltype, paste(prefix1, "cell_metadata.UMAPcluster.singleR_annot_cell_types.txt", sep="_"), 
            sep = "\t", quote = F, row.names = F)
print("done")
print("END OF SCRIPT")
