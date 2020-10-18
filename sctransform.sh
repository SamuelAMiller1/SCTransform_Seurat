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
ht29_ev <- PercentageFeatureSet(ht29_ev, pattern = "^MT-", col.name = "percent.mt"

ht29_lsd <- PercentageFeatureSet(ht29_lsd, pattern = "^MT-", col.name = "percent.mt"

##run sctransform
ht29_ev <- SCTransform(ht29_ev, vars.to.regress = "percent.mt", verbose = FALSE)

ht29_lsd <- SCTransform(ht29_lsd, vars.to.regress = "percent.mt", verbose = FALSE)


