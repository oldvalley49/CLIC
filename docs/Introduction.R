## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, eval=FALSE--------------------------------------------------------
# devtools::install_github("oldvalley49/CLIC")

## ----eval=FALSE---------------------------------------------------------------
# library(CLIC)
# library(Seurat)
# library(Signac)
# library(ggplot2)
# # replace with local path to data
# load("data/BMMC-s4d8.Rdata")
# 
# # RNA preprocessing
# rna <- NormalizeData(rna)
# out <- FindCLICFeatures(rna, species='human-signac-pearson')
# use_features <- out$use_features
# head(use_features)
# rna <- ScaleData(rna, features=use_features)
# # ATAC preprocessing
# atac <- RunTFIDF(atac)
# atac <- FindTopFeatures(atac, min.cutoff = "q0")
# atac <- RunSVD(atac)
# 
# DefaultAssay(atac) <- "ACTIVITY"
# atac <- NormalizeData(atac)
# atac <- ScaleData(atac, features = use_features)

## ----eval=FALSE---------------------------------------------------------------
# # identify anchors for cca
# transfer.anchors <- FindTransferAnchors(reference=rna, query=atac,
#                                         features=use_features,
#                                         reference.assay = "RNA",
#                                         query.assay = "ACTIVITY",
#                                         reduction = "cca")
# # coembed into shared latent space
# refdata <- GetAssayData(rna, assay = "RNA", layer = "data")
# imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]],
#                             dims = 2:30)
# atac[["RNA"]] <- imputation
# coembed <- merge(x = rna, y = atac, add.cell.ids = c("RNA", "ATAC"))
# coembed <- ScaleData(coembed, features = use_features, do.scale = FALSE)
# coembed <- RunPCA(coembed, features = use_features, verbose = FALSE)
# coembed <- RunUMAP(coembed, dims = 1:30)
# embedding_filtered <- Embeddings(coembed, reduction = "pca")
# p <- DimPlot(coembed, group.by = c("orig.ident"))

## ----include=FALSE, eval=FALSE------------------------------------------------
# ggsave("data/coembed_filtered.jpeg")

## ----echo=FALSE, out.width="500px"--------------------------------------------
knitr::include_graphics(system.file("extdata", "coembed_filtered.jpeg", package="CLIC"))

