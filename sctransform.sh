#################################
###Using sctransform in seurat###
#################################

##bash header
#!/bin/bash

## activate scRNAseq enviornment on miniconda
conda activate scRNAseq

##Navigate to output directory
cd /N/slate/millesaa/ht29_scRNA_seq/sctransform/

## load R
R

##Load data outputted from cell ranger and then create Seurat object
ht29_ev_data <- Read10X(data.dir = "/N/slate/millesaa/ht29_scRNAseq/ht29_ev/outs/filtered_feature_bc_matrix")

ht29_ev <- CreateSeuratObject(counts = ht29_ev_data)

ht29_lsd_data <- Read10X(data.dir = "/N/slate/millesaa/ht29_scRNAseq/ht29_lsd/outs/filtered_feature_bc_matrix")

ht29_lsd <- CreateSeuratObject(counts = ht29_lsd_data)

##store mitochondrial percentage in object meta data
ht29_ev <- PercentageFeatureSet(ht29_ev, pattern = "^MT-", col.name = "percent.mt")

ht29_lsd <- PercentageFeatureSet(ht29_lsd, pattern = "^MT-", col.name = "percent.mt")

##run sctransform
ht29_ev <- SCTransform(ht29_ev, vars.to.regress = "percent.mt", verbose = FALSE)

ht29_lsd <- SCTransform(ht29_lsd, vars.to.regress = "percent.mt", verbose = FALSE)

##visualization and clustering
ht29_ev <- RunPCA(ht29_ev, verbose = FALSE)
ht29_ev <- RunUMAP(ht29_ev, dims = 1:30, verbose = FALSE)

ht29_ev <- FindNeighbors(ht29_ev, dims = 1:30, verbose = FALSE)
ht29_ev <- FindClusters(ht29_ev, verbose = FALSE)
DimPlot(ht29_ev, label = TRUE) + NoLegend()

ht29_lsd <- RunPCA(ht29_lsd, verbose = FALSE) 
ht29_lsd <- RunUMAP(ht29_lsd, dims = 1:30, verbose = FALSE)

ht29_lsd <- FindNeighbors(ht29_lsd, dims = 1:30, verbose = FALSE)
ht29_lsd <- FindClusters(ht29_lsd, verbose = FALSE)
DimPlot(ht29_lsd, label = TRUE) + NoLegend()

##visualize canonical marker genes as violin plots
VlnPlot(ht29_ev, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"), 
    pt.size = 0.2, ncol = 4)

VlnPlot(ht29_lsd, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"), 
    pt.size = 0.2, ncol = 4)

##visualize canonical marker genes on the sctransform embedding
FeaturePlot(ht29_ev, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2, 
    ncol = 3)

FeaturePlot(ht29_lsd, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2, 
    ncol = 3)

