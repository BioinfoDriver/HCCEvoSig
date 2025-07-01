
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)

cellInfo <- read.csv2(file = '/data/GSE189903/GSE189903_Info.txt.gz', 
                      header = TRUE, sep = '\t')
cellInfo <- cellInfo %>% column_to_rownames(var = 'Cell')


expMatrix <- Read10X(data.dir = "/data/GSE189903")
multiHccSC <- CreateSeuratObject(counts = expMatrix, project = "multiHcc", min.cells = 10, min.features = 500)

# The percentage of reads that map to the mitochondrial genome
multiHccSC <- PercentageFeatureSet(object = multiHccSC, pattern = "^MT-", col.name = "perMito")

multiHccSC <- AddMetaData(object = multiHccSC, metadata = cellInfo)
multiHccSC <- subset(multiHccSC, perMito < 20)
multiHccSC <- subset(multiHccSC, nFeature_RNA < 6000)
multiHccSC <- subset(multiHccSC, nCount_RNA < 50000)
# 110817

# multiHccSC <- subset(multiHccSC, Type != '')


# Normalizing the data
multiHccSC <- NormalizeData(multiHccSC, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
multiHccSC <- FindVariableFeatures(multiHccSC, selection.method = "vst", nfeatures = 2000)

# Scaling the data
multiHccSC <- ScaleData(multiHccSC, features = rownames(multiHccSC))

# Perform linear dimensional reduction
multiHccSC <- RunPCA(multiHccSC, features = VariableFeatures(object = multiHccSC))

# Cluster the cells
multiHccSC <- FindNeighbors(object = multiHccSC, reduction = "pca", dims = 1:20)
multiHccSC <- FindClusters(multiHccSC, resolution = 0.8)

# Run non-linear dimensional reduction (UMAP/tSNE)
multiHccSC <- RunUMAP(multiHccSC, reduction = "pca", dims = 1:20)
multiHccSC <- RunTSNE(multiHccSC, reduction = "pca", dims = 1:20)



# Calculate module scores for feature expression programs
cellMarkers <- list(Tcells = c('CD2', 'CD3E', 'CD3D', 'CD3G'), 
                    Bcells = c('CD79A', 'SLAMF7', 'BLNK', 'FCRL5'), 
                    TECs = c('PECAM1', 'VWF', 'ENG', 'CDH5'), 
                    CAFs = c('COL1A2', 'FAP', 'PDPN', 'DCN', 'COL3A1', 'COL6A1'), 
                    TAMs = c('CD14', 'CD163', 'CD68', 'CSF1R'),  
                    EPIs = c('KRT8', 'KRT15', 'KRT16', 'KRT17', 'KRT18', 'KRT19', 'SFN', 'KRTCAP3', 'EPCAM'),
                    HEPs = c('ALB', 'TF', 'TTR'), CHOs = c('EPCAM', 'CD24', 'KRT19'))

multiHccSC <- AddModuleScore(object = multiHccSC, features = cellMarkers)

multiHccSC@meta.data <- multiHccSC@meta.data %>% 
  rename(Tcells = Cluster1, Bcells = Cluster2, TECs = Cluster3, 
         CAFs = Cluster4, TAMs = Cluster5, EPIs = Cluster6, HEPs = Cluster7, CHOs = Cluster8)



features <- c('Tcells', 'Bcells', 'TECs', 'CAFs', 'TAMs', 'EPIs', 'HEPs', 'CHOs')
pdf('/result/Section3/aver_lineage_specific_markers.pdf', height = 8, width = 16)
FeaturePlot(multiHccSC, features = features, reduction = 'umap', ncol = 4, raster = T, min.cutoff = "q1", max.cutoff = "q90")
FeaturePlot(multiHccSC, features = features, reduction = 'tsne', ncol = 4, raster = T, min.cutoff = "q1", max.cutoff = "q90")
dev.off()



# Cluster distribution
pdf('/result/Section3/view_cluster_distribution.pdf', height = 16, width = 16)
cowplot::plot_grid(ncol = 2, nrow = 2, 
  DimPlot(multiHccSC, reduction = 'umap', raster = T, label = T) + NoLegend(),
  DimPlot(multiHccSC, reduction = 'tsne', raster = T, label = T) + NoLegend(),
  DimPlot(multiHccSC, reduction = 'umap', raster = T, label = T, group.by = 'Type') + NoLegend(),
  DimPlot(multiHccSC, reduction = 'tsne', raster = T, label = T, group.by = 'Type') + NoLegend())
dev.off()



clusterIDs <- c(setNames(rep('Tcells', 13), c(0:8, 11, 18, 21, 24)), 
                setNames(rep('Bcells', 3), c(12, 17, 25)), 
                setNames(rep('TECs', 1), c(14)), 
                setNames(rep('CAFs', 1), c(13)), 
                setNames(rep('TAMs', 4), c(9, 10, 15, 16)),
                setNames(rep('EPIs', 4), c(19, 20, 22, 23)))                


# Assigning cell type identity to clusters
# Idents(multiHccSC) <- 'RNA_snn_res.0.8'
# multiHccSC <- SeuratObject::RenameIdents(multiHccSC, clusterIDs)
# multiHccSC$majorCellTypes <- Idents(multiHccSC)
# multiHccSC@meta.data <- multiHccSC@meta.data %>% rename(majorClusters = seurat_clusters)


multiHccSC@meta.data <- multiHccSC@meta.data %>% mutate(majorCellTypes = recode(seurat_clusters, !!!clusterIDs))
multiHccSC@meta.data <- multiHccSC@meta.data %>% rename(majorClusters = seurat_clusters)
Idents(multiHccSC) <- 'majorCellTypes'


###############################################Epis
Epis <- subset(multiHccSC, majorCellTypes == 'EPIs')
# 110817

# cycle
{
# Normalizing the data
Epis <- NormalizeData(Epis, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
Epis <- FindVariableFeatures(Epis, selection.method = "vst", nfeatures = 2000)

# Scaling the data
Epis <- ScaleData(Epis, features = rownames(Epis))

# Perform linear dimensional reduction
Epis <- RunPCA(Epis, features = VariableFeatures(object = Epis))

# Cluster the cells
Epis <- FindNeighbors(object = Epis, reduction = "pca", dims = 1:20)
Epis <- FindClusters(Epis, resolution = 0.2)

# Run non-linear dimensional reduction (UMAP/tSNE)
Epis <- RunUMAP(Epis, reduction = "pca", dims = 1:20)
Epis <- RunTSNE(Epis, reduction = "pca", dims = 1:20)
}


features <- c('Tcells', 'Bcells', 'TECs', 'CAFs', 'TAMs', 'EPIs', 'HEPs', 'CHOs')
pdf('/result/Section3/epis_aver_lineage_specific_markers.pdf', height = 8, width = 16)
FeaturePlot(Epis, features = features, reduction = 'umap', ncol = 4, raster = F, min.cutoff = "q1", max.cutoff = "q90")
FeaturePlot(Epis, features = features, reduction = 'tsne', ncol = 4, raster = F, min.cutoff = "q1", max.cutoff = "q90")
dev.off()


pdf('/result/Section3/epis_lineage_specific_markers.pdf', height = 8, width = 12)
FeaturePlot(Epis, features = c('KRT8', 'KRT18', 'KRT19', 'EPCAM', 'COL6A1', 'DCN'), 
            reduction = 'umap', ncol = 3, raster = F, min.cutoff = "q1", max.cutoff = "q90")
dev.off()



# Cluster distribution
pdf('/result/Section3/view_epis_cluster_distribution.pdf', height = 16, width = 16)
cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DimPlot(Epis, reduction = 'umap', raster = F, label = T) + NoLegend(),
                   DimPlot(Epis, reduction = 'tsne', raster = F, label = T) + NoLegend(),
                   DimPlot(Epis, reduction = 'umap', raster = F, label = T, group.by = 'Type') + NoLegend(),
                   DimPlot(Epis, reduction = 'tsne', raster = F, label = T, group.by = 'Type') + NoLegend())
dev.off()


# cycle
`%notin%` <- Negate(`%in%`)
Epis <- subset(Epis, RNA_snn_res.0.2 %notin% c(2, 7))
#



###############################################Tcells
Tcells <- subset(multiHccSC, majorCellTypes == 'Tcells')
# 93711

# cycle
{
  # Normalizing the data
  Tcells <- NormalizeData(Tcells, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identification of highly variable features (feature selection)
  Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000)
  
  # Scaling the data
  Tcells <- ScaleData(Tcells, features = rownames(Tcells))
  
  # Perform linear dimensional reduction
  Tcells <- RunPCA(Tcells, features = VariableFeatures(object = Tcells))
  
  # Cluster the cells
  Tcells <- FindNeighbors(object = Tcells, reduction = "pca", dims = 1:20)
  Tcells <- FindClusters(Tcells, resolution = 0.8)
  
  # Run non-linear dimensional reduction (UMAP/tSNE)
  Tcells <- RunUMAP(Tcells, reduction = "pca", dims = 1:20)
  Tcells <- RunTSNE(Tcells, reduction = "pca", dims = 1:20)
}


features <- c('Tcells', 'Bcells', 'TECs', 'CAFs', 'TAMs', 'EPIs', 'HEPs', 'CHOs')
pdf('/result/Section3/Tcells_aver_lineage_specific_markers.pdf', height = 8, width = 16)
FeaturePlot(Tcells, features = features, reduction = 'umap', ncol = 4, raster = F, min.cutoff = "q1", max.cutoff = "q90")
# FeaturePlot(Tcells, features = features, reduction = 'tsne', ncol = 4, raster = T, min.cutoff = "q1", max.cutoff = "q90")
dev.off()


pdf('/result/Section3/Tcells_lineage_specific_markers.pdf', height = 8, width = 8)
FeaturePlot(Tcells, features = c('CD2', 'CD3E', 'CD3D', 'CD3G'), reduction = 'umap', ncol = 2, raster = F, min.cutoff = "q1", max.cutoff = "q90")
dev.off()


# Cluster distribution
pdf('/result/Section3/view_Tcells_cluster_distribution.pdf', height = 16, width = 16)
cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DimPlot(Tcells, reduction = 'umap', raster = T, label = T) + NoLegend(),
                   DimPlot(Tcells, reduction = 'tsne', raster = T, label = T) + NoLegend(),
                   DimPlot(Tcells, reduction = 'umap', raster = T, label = T, group.by = 'Type') + NoLegend(),
                   DimPlot(Tcells, reduction = 'tsne', raster = T, label = T, group.by = 'Type') + NoLegend())
dev.off()


tcellIndex <- data.frame(TcellExp = colSums(GetAssayData(Tcells[["RNA"]], slot = "counts")[c('CD3D', 'CD3E', 'CD3G'), ] > 0) > 0)
tcellIndex <- merge(tcellIndex, Tcells[[]][, c('Type', 'RNA_snn_res.0.8')], by = 'row.names')


# cycle
`%notin%` <- Negate(`%in%`)
Tcells <- subset(Tcells, RNA_snn_res.0.8 %notin% c(9, 14, 15, 17, 20))
#



###############################################TAMs
TAMs <- subset(multiHccSC, majorCellTypes == 'TAMs')
# 8346

# cycle
{
  # Normalizing the data
  TAMs <- NormalizeData(TAMs, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identification of highly variable features (feature selection)
  TAMs <- FindVariableFeatures(TAMs, selection.method = "vst", nfeatures = 2000)
  
  # Scaling the data
  TAMs <- ScaleData(TAMs, features = rownames(TAMs))
  
  # Perform linear dimensional reduction
  TAMs <- RunPCA(TAMs, features = VariableFeatures(object = TAMs))
  
  # Cluster the cells
  TAMs <- FindNeighbors(object = TAMs, reduction = "pca", dims = 1:20)
  TAMs <- FindClusters(TAMs, resolution = 0.2)
  
  # Run non-linear dimensional reduction (UMAP/tSNE)
  TAMs <- RunUMAP(TAMs, reduction = "pca", dims = 1:20)
  TAMs <- RunTSNE(TAMs, reduction = "pca", dims = 1:20)
}


features <- c('Tcells', 'Bcells', 'TECs', 'CAFs', 'TAMs', 'EPIs', 'HEPs', 'CHOs')
pdf('/result/Section3/tams_aver_lineage_specific_markers.pdf', height = 8, width = 16)
FeaturePlot(TAMs, features = features, reduction = 'umap', ncol = 4, raster = F, min.cutoff = "q1", max.cutoff = "q90")
FeaturePlot(TAMs, features = features, reduction = 'tsne', ncol = 4, raster = F, min.cutoff = "q1", max.cutoff = "q90")
dev.off()


pdf('/result/Section3/TAMs_lineage_specific_markers.pdf', height = 8, width = 12)
FeaturePlot(TAMs, features = c('CD14', 'CD163', 'CD68', 'CSF1R', 'SLAMF7', 'CD19'), 
            reduction = 'umap', ncol = 3, raster = F, min.cutoff = "q1", max.cutoff = "q90")
dev.off()



# Cluster distribution
pdf('/result/Section3/view_tams_cluster_distribution.pdf', height = 16, width = 16)
cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DimPlot(TAMs, reduction = 'umap', raster = F, label = T) + NoLegend(),
                   DimPlot(TAMs, reduction = 'tsne', raster = F, label = T) + NoLegend(),
                   DimPlot(TAMs, reduction = 'umap', raster = F, label = T, group.by = 'Type') + NoLegend(),
                   DimPlot(TAMs, reduction = 'tsne', raster = F, label = T, group.by = 'Type') + NoLegend())
dev.off()


# cycle
`%notin%` <- Negate(`%in%`)
TAMs <- subset(TAMs, RNA_snn_res.0.2 %notin% c(5))
TAMs <- subset(TAMs, RNA_snn_res.0.2 %notin% c(6))
#



###############################################fMultiHccSC
fMultiHccSC <- subset(multiHccSC, cells = c(colnames(Epis), colnames(TAMs), colnames(Tcells),
                                            rownames(subset(multiHccSC[[]], majorCellTypes %in% c('Bcells', 'CAFs', 'TECs')))))

# Normalizing the data
fMultiHccSC <- NormalizeData(fMultiHccSC, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
fMultiHccSC <- FindVariableFeatures(fMultiHccSC, selection.method = "vst", nfeatures = 2000)

# Scaling the data
fMultiHccSC <- ScaleData(fMultiHccSC, features = rownames(fMultiHccSC))

# Perform linear dimensional reduction
fMultiHccSC <- RunPCA(fMultiHccSC, features = VariableFeatures(object = fMultiHccSC))

# Cluster the cells
fMultiHccSC <- FindNeighbors(object = fMultiHccSC, reduction = "pca", dims = 1:30)
fMultiHccSC <- FindClusters(fMultiHccSC, resolution = 0.8)

# Run non-linear dimensional reduction (UMAP/tSNE)
fMultiHccSC <- RunUMAP(fMultiHccSC, reduction = "pca", dims = 1:30)
fMultiHccSC <- RunTSNE(fMultiHccSC, reduction = "pca", dims = 1:30)


pdf('/result/Section3/fMultiHccSC_ElbowPlot.pdf')
ElbowPlot(fMultiHccSC, ndims = 50)

dev.off()



features <- c('Tcells', 'Bcells', 'TECs', 'CAFs', 'TAMs', 'EPIs', 'HEPs', 'CHOs')
pdf('/result/Section3/fMultiHccSC_aver_lineage_specific_markers.pdf', height = 8, width = 16)
FeaturePlot(fMultiHccSC, features = features, reduction = 'umap', ncol = 4, raster = T, min.cutoff = "q1", max.cutoff = "q90", cols = c("lightgrey", "#ff0000"))
FeaturePlot(fMultiHccSC, features = features, reduction = 'tsne', ncol = 4, raster = T, min.cutoff = "q1", max.cutoff = "q90", cols = c("lightgrey", "#ff0000"))
dev.off()



# pdf('/result/Section3/fMultiHccSC_lineage_specific_markers.pdf', height = 12, width = 8)
# plots <- FeaturePlot(fMultiHccSC, features = c('CD2', 'CD3E', 'CD3D', 'CD3G', 'CD79A', 'SLAMF7', 'BLNK', 'FCRL5', 'PECAM1', 'VWF', 'ENG', 'CDH5', 
#                                'DCN', 'COL3A1', 'COL6A1', 'FAP', 'CD14', 'CD163', 'CD68', 'CSF1R', 'EPCAM', 'KRT8', 'KRT18', 'KRT19'), 
#             reduction = 'umap', ncol = 4, raster = T, min.cutoff = "q1", max.cutoff = "q90", cols = c("lightgrey", "#ff0000"), combine = FALSE)
# 
# plots <- lapply(plots, function(p) p + NoLegend())
# cowplot::plot_grid(plotlist = plots, ncol = 4)
# 
# dev.off()


pdf('/result/Section3/fMultiHccSC_lineage_specific_markers.pdf', height = 8, width = 12)
plots <- FeaturePlot(fMultiHccSC, features = c('CD3E', 'CD79A', 'PECAM1', 'DCN', 'CD68', 'KRT8'), 
                     reduction = 'umap', ncol = 3, raster = T, min.cutoff = "q1", max.cutoff = "q90", 
                     cols = c("lightgrey", "#ff0000"), combine = FALSE, raster.dpi = c(256, 256), order = T)

# plots <- lapply(plots, function(p) p + NoLegend())
cowplot::plot_grid(plotlist = plots, ncol = 3)

dev.off()


# Cluster distribution
pdf('/result/Section3/view_fMultiHccSC_cluster_distribution.pdf', height = 8, width = 16)
cowplot::plot_grid(ncol = 2, nrow = 1, 
                   DimPlot(fMultiHccSC, reduction = 'umap', raster = T, label = T) + NoLegend(),
                   DimPlot(fMultiHccSC, reduction = 'tsne', raster = T, label = T) + NoLegend())
dev.off()


# saveRDS(fMultiHccSC, file = '/data/GSE189903/fMultiHccSC.rds')


pdf('/result/Section3/fMultiHccSC_cellMarkers_exp.pdf', height = 24, width = 32)

cellMarkers <- c('CD2', 'CD3E', 'CD3D', 'CD3G', 'CD79A', 'SLAMF7', 'BLNK', 'FCRL5', 'PECAM1', 'VWF', 'ENG', 'CDH5', 
                 'COL1A2', 'FAP', 'PDPN', 'DCN', 'COL3A1', 'COL6A1', 'CD14', 'CD163', 'CD68', 'CSF1R',  
                 'KRT8', 'KRT15', 'KRT16', 'KRT17', 'KRT18', 'KRT19', 'SFN', 'KRTCAP3', 'EPCAM',
                 'ALB', 'TF', 'TTR', 'CD24')

DotPlot(fMultiHccSC, features = cellMarkers, cols = c("lightgrey", "#ff0000"))
VlnPlot(fMultiHccSC, features = cellMarkers, pt.size = 0, ncol = 5)

dev.off()



# Assigning cell type identity to clusters

clusterIDs <- c(setNames(rep('Tcells', 18), c(0:7, 10, 12, 16, 18, 20, 22, 24, 26, 28, 32)), 
                setNames(rep('Bcells', 3), c(11, 17, 29)), 
                setNames(rep('TECs', 2), c(14, 30)), 
                setNames(rep('CAFs', 1), c(13)), 
                setNames(rep('TAMs', 4), c(8, 9, 15, 21)),
                setNames(rep('EPIs', 5), c(19, 23, 25, 27, 31)))                

fMultiHccSC@meta.data <- fMultiHccSC@meta.data %>% mutate(cellTypes = recode(seurat_clusters, !!!clusterIDs))
fMultiHccSC@meta.data <- fMultiHccSC@meta.data %>% rename(clusters = seurat_clusters)
Idents(fMultiHccSC) <- 'cellTypes'



# Cluster distribution
pdf('/result/Section3/view_fMultiHccSC_cellType_distribution.pdf', height = 8, width = 8)
cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DimPlot(fMultiHccSC, reduction = 'umap', raster = T, label = T, raster.dpi = c(256, 256)) + NoLegend(),
                   DimPlot(fMultiHccSC, reduction = 'tsne', raster = T, label = T, raster.dpi = c(256, 256)) + NoLegend())
dev.off()



pdf('/result/Section3/fMultiHccSC_cellType_cellMarkers_exp.pdf', height = 4, width = 8)

DotPlot(fMultiHccSC, features = c('CD3E', 'CD3D', 'CD3G', 'CD79A', 'SLAMF7', 'FCRL5', 'PECAM1', 'VWF', 'CDH5', 
                                  'COL1A2', 'DCN', 'COL3A1', 'CD14', 'CD68', 'CSF1R', 'KRT8', 'KRT18', 'EPCAM'), cols = c("lightgrey", "#ff0000")) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()




pdf('/result/Section3/fMultiHccSC_cellType_com.pdf', height = 8, width = 8)

comdata <- filterfMultiHccSC[[]][, c('Type', 'cellTypes')] %>% 
  mutate(Type = ifelse(cellTypes == 'EPIs' & Type == 'unclassified', 'non-malignant cell', Type)) %>% 
  subset.data.frame(Type != 'unclassified') %>% 
  group_by(cellTypes) %>% count(Type, name = 'Cells') %>% as.data.frame()

ggplot(data = comdata, aes(x = cellTypes, y = Cells, fill = Type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = NULL, y = "Percentage", fill = "Cell types") + theme_minimal()

dev.off()




sigs <- c('ADH4', 'CDC20', 'CFHR3', 'CYP2C9', 'RAMP3', 'RDH16', 'SERPINE1', 'SLC16A11', 'SPP1', 'SPP2', 'TRNP1')
colors<- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid", "orange", "gold", "gray", "#ff0000", "blue")
# colors<- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid")
sams <- c('1CT1', '1CT2', '1CT3', '1HT1', '1HT2', '1HT3', '2CT1', '2CT2', '2CT3', 
          '2HT1', '2HT2', '2HT3', '3CT1', '3CT2', '3HT1', '3HT2', '3HT3', '4HT1', '4HT2', '4HT3')


pdf('/result/Section3/fMultiHccSC_cellType_sigs_exp.pdf', height = 16, width = 16)

# RidgePlot(fMultiHccSC, features = sigs, ncol = 2)
cowplot::plot_grid(ncol = 2, nrow = 3, rel_heights = c(1, 2, 1), 
                   DotPlot(subset(fMultiHccSC, Sample %in% sams), features = sigs, cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(fMultiHccSC, features = sigs, cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   VlnPlot(subset(fMultiHccSC, Sample %in% sams), features = sigs, pt.size = 0, stack = TRUE, cols = colors, flip = T) + NoLegend(), 
                   VlnPlot(fMultiHccSC, features = sigs, pt.size = 0, stack = TRUE, cols = colors, flip = T) + NoLegend())

plots <- FeaturePlot(subset(fMultiHccSC, Sample %in% sams), features = sigs, 
                     reduction = 'umap', ncol = 3, raster = T, min.cutoff = "q5", max.cutoff = "q90", 
                     cols = c("lightgrey", "#ff0000"), combine = FALSE, raster.dpi = c(256, 256), order = T)

# plots <- lapply(plots, function(p) p + NoLegend())
cowplot::plot_grid(plotlist = plots, ncol = 3)


plots <- FeaturePlot(fMultiHccSC, features = sigs, 
                     reduction = 'umap', ncol = 3, raster = T, min.cutoff = "q5", max.cutoff = "q90", 
                     cols = c("lightgrey", "#ff0000"), combine = FALSE, raster.dpi = c(256, 256), order = T)

# plots <- lapply(plots, function(p) p + NoLegend())
cowplot::plot_grid(plotlist = plots, ncol = 3)

dev.off()




##################Analysis, a minimum of 10 cells
nCells <- fMultiHccSC[[]][, c('Sample', 'cellTypes')] %>% group_by(Sample, cellTypes) %>% count(name = 'nCells') %>% subset(nCells >=10)
nCells <- merge(fMultiHccSC[[]] %>% tibble::rownames_to_column(var = "cellNames"), nCells, by = c('Sample', 'cellTypes'), all = FALSE)
filterfMultiHccSC <- subset(fMultiHccSC, cells = nCells$cellNames)


# Diff expression
filterfMultiHccSC$sampleType <- factor(substr(filterfMultiHccSC$Sample, 3, 3), levels = c('N', 'B', 'T'))
filterfMultiHccSC$diseaseType <- factor(substr(filterfMultiHccSC$Sample, 2, 2), levels = c('C', 'H'))



pdf('/result/Section3/fMultiHccSC_cellType_sigs_HCCexp.pdf', height = 16, width = 16)

cowplot::plot_grid(ncol = 2, nrow = 3, rel_heights = c(1, 2, 1), 
                   DotPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), features = sigs, cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, diseaseType == 'H'), features = sigs, cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   VlnPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), features = sigs, pt.size = 0, stack = TRUE, cols = colors, flip = T) + NoLegend(), 
                   VlnPlot(subset(filterfMultiHccSC, diseaseType == 'H'), features = sigs, pt.size = 0, stack = TRUE, cols = colors, flip = T) + NoLegend())

plots <- FeaturePlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), features = sigs, 
                     reduction = 'umap', ncol = 3, raster = T, min.cutoff = "q5", max.cutoff = "q90", 
                     cols = c("lightgrey", "#ff0000"), combine = FALSE, raster.dpi = c(256, 256), order = T)

cowplot::plot_grid(plotlist = plots, ncol = 3)

plots <- FeaturePlot(subset(filterfMultiHccSC, diseaseType == 'H'), features = sigs, 
                     reduction = 'umap', ncol = 3, raster = T, min.cutoff = "q5", max.cutoff = "q90", 
                     cols = c("lightgrey", "#ff0000"), combine = FALSE, raster.dpi = c(256, 256), order = T)

cowplot::plot_grid(plotlist = plots, ncol = 3)

dev.off()



pdf('/result/Section3/fMultiHccSC_cellType_sigs_multiregion_exp.pdf', height = 16, width = 16)

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DotPlot(subset(filterfMultiHccSC, Sample %in% sams), idents = 'EPIs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, Sample %in% sams), idents = 'TECs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, Sample %in% sams), idents = 'CAFs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, Sample %in% sams), idents = 'TAMs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)))

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DotPlot(filterfMultiHccSC, idents = 'EPIs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(filterfMultiHccSC, idents = 'TECs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(filterfMultiHccSC, idents = 'CAFs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(filterfMultiHccSC, idents = 'TAMs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)))

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   VlnPlot(subset(filterfMultiHccSC, Sample %in% sams), features = sigs, idents = 'EPIs', group.by = 'Sample', pt.size = 0, stack = TRUE, flip = T) + NoLegend(), 
                   VlnPlot(subset(filterfMultiHccSC, Sample %in% sams), features = c('RAMP3', 'SERPINE1', 'PECAM1', 'VWF', 'ENG', 'CDH5', 'PLVAP', 'CD34', 'CLDN5', 'KDR', 'FLT1'), idents = 'TECs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(subset(filterfMultiHccSC, Sample %in% sams), features = c('SERPINE1','TRNP1', 'COL1A2', 'FAP', 'PDPN', 'DCN', 'COL3A1', 'COL6A1', 'PDGFRA', 'THY1', 'LUM'), idents = 'CAFs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(subset(filterfMultiHccSC, Sample %in% sams), features = c('SPP1', 'CDC20', 'CD14', 'CD163', 'CD68', 'CSF1R', 'AIF1', 'ADGRE1', 'MRC1', 'C1QB', 'APOE'), idents = 'TAMs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend())

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   VlnPlot(filterfMultiHccSC, features = sigs, idents = 'EPIs', group.by = 'Sample', pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(filterfMultiHccSC, features = c('RAMP3', 'SERPINE1', 'PECAM1', 'VWF', 'ENG', 'CDH5', 'PLVAP', 'CD34', 'CLDN5', 'KDR', 'FLT1'), idents = 'TECs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(filterfMultiHccSC, features = c('SERPINE1','TRNP1', 'COL1A2', 'FAP', 'PDPN', 'DCN', 'COL3A1', 'COL6A1', 'PDGFRA', 'THY1', 'LUM'), idents = 'CAFs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(filterfMultiHccSC, features = c('SPP1', 'CDC20', 'CD14', 'CD163', 'CD68', 'CSF1R', 'AIF1', 'ADGRE1', 'MRC1', 'C1QB', 'APOE'), idents = 'TAMs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend())

dev.off()
                   




pdf('/result/Section3/fMultiHccSC_cellType_sigs_HCCmultiregion_exp.pdf', height = 16, width = 16)

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DotPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), idents = 'EPIs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), idents = 'TECs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), idents = 'CAFs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), idents = 'TAMs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)))

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DotPlot(subset(filterfMultiHccSC, diseaseType == 'H'), idents = 'EPIs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, diseaseType == 'H'), idents = 'TECs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, diseaseType == 'H'), idents = 'CAFs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   DotPlot(subset(filterfMultiHccSC, diseaseType == 'H'), idents = 'TAMs', features = sigs, group.by = 'Sample', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)))

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   VlnPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), features = sigs, idents = 'EPIs', group.by = 'Sample', pt.size = 0, stack = TRUE, flip = T) + NoLegend(), 
                   VlnPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), features = c('RAMP3', 'SERPINE1', 'PECAM1', 'VWF', 'ENG', 'CDH5', 'PLVAP', 'CD34', 'CLDN5', 'KDR', 'FLT1'), idents = 'TECs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), features = c('SERPINE1','TRNP1', 'COL1A2', 'FAP', 'PDPN', 'DCN', 'COL3A1', 'COL6A1', 'PDGFRA', 'THY1', 'LUM'), idents = 'CAFs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(subset(filterfMultiHccSC, Sample %in% sams & diseaseType == 'H'), features = c('SPP1', 'CDC20', 'CD14', 'CD163', 'CD68', 'CSF1R', 'AIF1', 'ADGRE1', 'MRC1', 'C1QB', 'APOE'), idents = 'TAMs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend())

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   VlnPlot(subset(filterfMultiHccSC, diseaseType == 'H'), features = sigs, idents = 'EPIs', group.by = 'Sample', pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(subset(filterfMultiHccSC, diseaseType == 'H'), features = c('RAMP3', 'SERPINE1', 'PECAM1', 'VWF', 'ENG', 'CDH5', 'PLVAP', 'CD34', 'CLDN5', 'KDR', 'FLT1'), idents = 'TECs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(subset(filterfMultiHccSC, diseaseType == 'H'), features = c('SERPINE1','TRNP1', 'COL1A2', 'FAP', 'PDPN', 'DCN', 'COL3A1', 'COL6A1', 'PDGFRA', 'THY1', 'LUM'), idents = 'CAFs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend(),
                   VlnPlot(subset(filterfMultiHccSC, diseaseType == 'H'), features = c('SPP1', 'CDC20', 'CD14', 'CD163', 'CD68', 'CSF1R', 'AIF1', 'ADGRE1', 'MRC1', 'C1QB', 'APOE'), idents = 'TAMs', group.by = 'Sample',pt.size = 0, stack = TRUE, flip = T) + NoLegend())

dev.off()





##########Epithelial
pdf('/result/Section3/fMultiHccSC_sigs_epi_diffexp.pdf', height = 8, width = 8)

cowplot::plot_grid(ncol = 2, nrow = 2, 
  DoHeatmap(object = subset(filterfMultiHccSC, cellTypes == 'EPIs'), features = sigs, group.by = "sampleType", raster = T),
  
  DotPlot(filterfMultiHccSC, idents = 'EPIs', features = sigs, group.by = 'sampleType', cols = c("lightgrey", "#ff0000")) + 
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),

  VlnPlot(filterfMultiHccSC, features = sigs, idents = 'EPIs', group.by = 'sampleType', pt.size = 0, stack = TRUE, flip = T) + NoLegend())

VlnPlot(filterfMultiHccSC, features = sigs, idents = 'EPIs', group.by = 'sampleType', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()


pdf('/result/Section3/fMultiHccSC_sigs_HCCepi_diffexp.pdf', height = 8, width = 8)

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DoHeatmap(object = subset(filterfMultiHccSC, cellTypes == 'EPIs' & diseaseType == 'H'), features = sigs, group.by = "sampleType", raster = T),
                   
                   DotPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), idents = 'EPIs', features = sigs, group.by = 'sampleType', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   
                   VlnPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), features = sigs, idents = 'EPIs', group.by = 'sampleType', pt.size = 0, stack = TRUE, flip = T) + NoLegend())

VlnPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), features = sigs, idents = 'EPIs', group.by = 'sampleType', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()



##########Endothelial
pdf('/result/Section3/fMultiHccSC_sigs_tec_diffexp.pdf', height = 8, width = 8)

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DoHeatmap(object = subset(filterfMultiHccSC, cellTypes == 'TECs'), features = sigs, group.by = "sampleType", raster = T),
                   
                   DotPlot(filterfMultiHccSC, idents = 'TECs', features = sigs, group.by = 'sampleType', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   
                   VlnPlot(filterfMultiHccSC, features = sigs, idents = 'TECs', group.by = 'sampleType', pt.size = 0, stack = TRUE, flip = T) + NoLegend())

VlnPlot(filterfMultiHccSC, features = sigs, idents = 'TECs', group.by = 'sampleType', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()


pdf('/result/Section3/fMultiHccSC_sigs_HCCtec_diffexp.pdf', height = 8, width = 8)
tecSigs<- c(setdiff(sigs, 'CFHR3'), 'CDH5')

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DoHeatmap(object = subset(filterfMultiHccSC, cellTypes == 'TECs' & diseaseType == 'H'), features = tecSigs, group.by = "sampleType", raster = T),
                   
                   DotPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), idents = 'TECs', features = tecSigs, group.by = 'sampleType', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   
                   VlnPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), features = tecSigs, idents = 'TECs', group.by = 'sampleType', pt.size = 0, stack = TRUE, flip = T) + NoLegend())

VlnPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), features = tecSigs, idents = 'TECs', group.by = 'sampleType', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()



##########Fibroblasts
pdf('/result/Section3/fMultiHccSC_sigs_caf_diffexp.pdf', height = 8, width = 8)

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DoHeatmap(object = subset(filterfMultiHccSC, cellTypes == 'CAFs'), features = sigs, group.by = "sampleType", raster = T),
                   
                   DotPlot(filterfMultiHccSC, idents = 'CAFs', features = sigs, group.by = 'sampleType', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   
                   VlnPlot(filterfMultiHccSC, features = sigs, idents = 'CAFs', group.by = 'sampleType', pt.size = 0, stack = TRUE, flip = T) + NoLegend())

VlnPlot(filterfMultiHccSC, features = sigs, idents = 'CAFs', group.by = 'sampleType', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()


pdf('/result/Section3/fMultiHccSC_sigs_HCCcaf_diffexp.pdf', height = 8, width = 8)
cafSigs<- c(setdiff(sigs, 'RDH16'), 'DCN')

cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DoHeatmap(object = subset(filterfMultiHccSC, cellTypes == 'CAFs' & diseaseType == 'H'), features = cafSigs, group.by = "sampleType", raster = T),
                   
                   DotPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), idents = 'CAFs', features = cafSigs, group.by = 'sampleType', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   
                   VlnPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), features = cafSigs, idents = 'CAFs', group.by = 'sampleType', pt.size = 0, stack = TRUE, flip = T) + NoLegend())

VlnPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), features = cafSigs, idents = 'CAFs', group.by = 'sampleType', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()




##########Macrophages
pdf('/result/Section3/fMultiHccSC_sigs_tam_diffexp.pdf', height = 8, width = 8)

tamSigs<- c('SPP1', 'CDC20', 'CD14', 'CD163', 'CD68', 'CSF1R', 'AIF1', 'ADGRE1', 'MRC1', 'C1QB', 'APOE')
cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DoHeatmap(object = subset(filterfMultiHccSC, cellTypes == 'TAMs'), features = tamSigs, group.by = "sampleType", raster = T),
                   
                   DotPlot(filterfMultiHccSC, idents = 'TAMs', features = tamSigs, group.by = 'sampleType', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   
                   VlnPlot(filterfMultiHccSC, features = tamSigs, idents = 'TAMs', group.by = 'sampleType', pt.size = 0, stack = TRUE, flip = T) + NoLegend())

VlnPlot(filterfMultiHccSC, features = tamSigs, idents = 'TAMs', group.by = 'sampleType', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()



pdf('/result/Section3/fMultiHccSC_sigs_HCCtam_diffexp.pdf', height = 8, width = 8)

tamSigs<- c('SPP1', 'CDC20', 'CD14', 'CD163', 'CD68', 'CSF1R', 'AIF1', 'ADGRE1', 'MRC1', 'C1QB', 'APOE')
cowplot::plot_grid(ncol = 2, nrow = 2, 
                   DoHeatmap(object = subset(filterfMultiHccSC, cellTypes == 'TAMs' & diseaseType == 'H'), features = tamSigs, group.by = "sampleType", raster = T),
                   
                   DotPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), idents = 'TAMs', features = tamSigs, group.by = 'sampleType', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   
                   VlnPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), features = tamSigs, idents = 'TAMs', group.by = 'sampleType', pt.size = 0, stack = TRUE, flip = T) + NoLegend())

VlnPlot(object = subset(filterfMultiHccSC, diseaseType == 'H'), features = tamSigs, idents = 'TAMs', group.by = 'sampleType', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()



##########Malignant cell
pdf('/result/Section3/fMultiHccSC_sigs_Malignantepi_diffexp.pdf', height = 8, width = 8)

cowplot::plot_grid(ncol = 1, nrow = 2, 
                   DotPlot(object = subset(filterfMultiHccSC, cellTypes == 'EPIs' & Type %in% c('Malignant cell', 'unclassified')), features = sigs, group.by = 'Type', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   
                   VlnPlot(object =  subset(filterfMultiHccSC, cellTypes == 'EPIs' & Type %in% c('Malignant cell', 'unclassified')), features = sigs, group.by = 'Type', pt.size = 0, stack = TRUE, flip = T) + NoLegend())


VlnPlot(object =  subset(filterfMultiHccSC, cellTypes == 'EPIs' & Type %in% c('Malignant cell', 'unclassified')), features = sigs, group.by = 'Type', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()



pdf('/result/Section3/fMultiHccSC_sigs_HCCMalignantepi_diffexp.pdf', height = 8, width = 8)

cowplot::plot_grid(ncol = 1, nrow = 2, 
                   DotPlot(object = subset(filterfMultiHccSC, cellTypes == 'EPIs' & Type %in% c('Malignant cell', 'unclassified') & diseaseType == 'H'), features = sigs, group.by = 'Type', cols = c("lightgrey", "#ff0000")) + 
                     theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
                   
                   VlnPlot(object =  subset(filterfMultiHccSC, cellTypes == 'EPIs' & Type %in% c('Malignant cell', 'unclassified') & diseaseType == 'H'), features = sigs, group.by = 'Type', pt.size = 0, stack = TRUE, flip = T) + NoLegend())


VlnPlot(object =  subset(filterfMultiHccSC, cellTypes == 'EPIs' & Type %in% c('Malignant cell', 'unclassified') & diseaseType == 'H'), features = sigs, group.by = 'Type', pt.size = 0, ncol = 4, flip = T) + NoLegend()
dev.off()



# differentially expressed genes
FindMarkers(object = filterfMultiHccSC, ident.1 = 'T', ident.2 = 'N', group.by = 'sampleType', min.pct = 0.05,
            subset.ident = 'EPIs', features = sigs, logfc.threshold = 0.05)

#               p_val avg_log2FC pct.1 pct.2    p_val_adj
# RDH16  2.871904e-82  1.1876322 0.359 0.027 6.572065e-78
# SPP2   2.433422e-33  0.9033112 0.435 0.117 5.568643e-29
# ADH4   1.025250e-31  1.1366594 0.555 0.193 2.346182e-27
# CYP2C9 2.125856e-23  1.1670595 0.498 0.214 4.864808e-19
# CFHR3  4.814007e-19  0.4568641 0.144 0.024 1.101637e-14
# TRNP1  9.631155e-06 -0.2798780 0.029 0.136 2.203994e-01
# SPP1   1.328421e-01 -1.6022971 0.196 0.225 1.000000e+00

FindMarkers(object = filterfMultiHccSC, ident.1 = 'T', ident.2 = 'N', group.by = 'sampleType', min.pct = 0.05,
            subset.ident = 'TECs', features = sigs, logfc.threshold = 0.05)

#                 p_val avg_log2FC pct.1 pct.2   p_val_adj
# RAMP3    1.638716e-07  0.5980569 0.643 0.478 0.003750038
# SERPINE1 8.433128e-03 -0.3874585 0.084 0.145 1.000000000

FindMarkers(object = filterfMultiHccSC, ident.1 = 'T', ident.2 = 'N', group.by = 'sampleType', min.pct = 0.05,
            subset.ident = 'CAFs', features = sigs, logfc.threshold = 0.05)
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# SERPINE1 0.004970742 -0.4576838 0.159 0.332         1

FindMarkers(object = filterfMultiHccSC, ident.1 = 'T', ident.2 = 'N', group.by = 'sampleType', min.pct = 0.05,
            subset.ident = 'TAMs', features = sigs, logfc.threshold = 0.05)
#              p_val avg_log2FC pct.1 pct.2     p_val_adj
# SPP1 5.230785e-191  -4.376879 0.008  0.41 1.197013e-186




#################HCC
# differentially expressed genes
FindMarkers(object = subset(filterfMultiHccSC, diseaseType == 'H'), ident.1 = 'T', ident.2 = 'N', group.by = 'sampleType', min.pct = 0.05,
            subset.ident = 'EPIs', features = sigs, logfc.threshold = 0.05)
#                 p_val avg_log2FC pct.1 pct.2    p_val_adj
# RDH16    1.995130e-67  1.3570974 0.406 0.025 4.565654e-63
# ADH4     8.192930e-25  1.0833961 0.650 0.234 1.874870e-20
# SPP2     2.175514e-24  0.9295152 0.497 0.142 4.978446e-20
# CYP2C9   2.793327e-18  1.1375551 0.566 0.255 6.392250e-14
# CFHR3    3.233416e-10  0.4849240 0.133 0.027 7.399348e-06
# SERPINE1 7.741928e-05  0.4203521 0.210 0.101 1.000000e+00
# CDC20    3.301619e-04 -0.5410407 0.021 0.122 1.000000e+00
# SPP1     3.343054e-02 -1.9311236 0.182 0.243 1.000000e+00


FindMarkers(object = subset(filterfMultiHccSC, diseaseType == 'H'), ident.1 = 'T', ident.2 = 'N', group.by = 'sampleType', min.pct = 0.05,
            subset.ident = 'TECs', features = sigs, logfc.threshold = 0.05)
#              p_val avg_log2FC pct.1 pct.2   p_val_adj
# RAMP3 9.393039e-08  0.8525513 0.672 0.483 0.002149503
# SPP1  5.292485e-04 -1.0016861 0.022 0.106 1.000000000


FindMarkers(object = subset(filterfMultiHccSC, diseaseType == 'H'), ident.1 = 'T', ident.2 = 'N', group.by = 'sampleType', min.pct = 0.05,
            subset.ident = 'CAFs', features = sigs, logfc.threshold = 0.05)
#            p_val avg_log2FC pct.1 pct.2 p_val_adj
# SERPINE1 0.73273  0.6379561 0.132 0.163         1


FindMarkers(object = subset(filterfMultiHccSC, diseaseType == 'H'), ident.1 = 'T', ident.2 = 'N', group.by = 'sampleType', min.pct = 0.05,
            subset.ident = 'TAMs', features = sigs, logfc.threshold = 0.05)
#              p_val avg_log2FC pct.1 pct.2     p_val_adj
# SPP1 2.157061e-237  -4.768374 0.005 0.508 4.936218e-233




FindMarkers(object = subset(filterfMultiHccSC, Type %in% c('Malignant cell', 'unclassified')), 
            ident.1 = 'Malignant cell', ident.2 = 'unclassified', group.by = 'Type', 
            min.pct = 0.05, subset.ident = 'EPIs', features = sigs, logfc.threshold = 0.05)
#                  p_val avg_log2FC pct.1 pct.2     p_val_adj
# ADH4     1.861608e-119 -2.2260158 0.065 0.492 4.260103e-115
# SPP2      5.255030e-93 -1.4143742 0.027 0.364  1.202561e-88
# CYP2C9    6.303147e-46 -1.4192512 0.162 0.401  1.442412e-41
# RDH16     3.876800e-40 -0.5314386 0.008 0.160  8.871669e-36
# SPP1      1.968978e-31  2.1294626 0.297 0.101  4.505809e-27
# SLC16A11  5.923490e-19 -0.7663067 0.117 0.248  1.355531e-14
# CDC20     9.002637e-15  0.3935487 0.128 0.034  2.060163e-10

FindMarkers(object = subset(filterfMultiHccSC, Type %in% c('Malignant cell', 'unclassified') & diseaseType == 'H'), 
            ident.1 = 'Malignant cell', ident.2 = 'unclassified', group.by = 'Type', 
            min.pct = 0.05, subset.ident = 'EPIs', features = sigs, logfc.threshold = 0.05)
#                 p_val avg_log2FC pct.1 pct.2    p_val_adj
# ADH4     1.952868e-88 -2.2173978 0.098 0.570 4.468943e-84
# SPP2     2.751803e-75 -1.5938792 0.039 0.440 6.297226e-71
# CYP2C9   3.546556e-43 -1.5940597 0.177 0.473 8.115939e-39
# RDH16    3.265336e-32 -0.6398384 0.006 0.184 7.472396e-28
# SPP1     2.144901e-25  2.1879701 0.337 0.111 4.908391e-21
# CDC20    9.784898e-15  0.4938061 0.162 0.034 2.239176e-10
# SLC16A11 8.833964e-11 -0.7824335 0.175 0.292 2.021564e-06


# saveRDS(filterfMultiHccSC, file = '/data/GSE189903/filterfMultiHccSC.rds')




