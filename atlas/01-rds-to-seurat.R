library(Seurat)
library(SeuratDisk)
atlas <- readRDS("scAtlas.rds")
#because SeuratDisk won't export counts to anndata.raw.X otherwise
dt <- Seurat::DietSeurat(atlas)
SeuratDisk::SaveH5Seurat(atlas, filename = "scAtlas.rds.h5Seurat", verbose = TRUE, overwrite = TRUE)
