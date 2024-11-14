install.packages("Seurat")
install.packages("SeuratObject")
library(Seurat)
#1:30, 1.5
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

library(SeuratData)
library(SeuratWrappers)
install.packages("harmony")

library("harmony")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("satijalab/seurat-wrappers")

BiocManager::install("batchelor")

library("batchelor")
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library("harmony")
library("batchelor")
library(ggplot2)
library(patchwork)
library(gridExtra)
library(grid)

options(future.globals.maxSize = 1e9) 
setwd("C:/Lab_Xu/Data/xu43_raw")

setwd("C:/Lab_Xu/Data/xu48_raw")

setwd("C:/Lab_Xu/Data/xu53_raw")

setwd("C:/Lab_Xu/Data/xu55_raw")

Maize.xu43FM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu43_raw/raw_feature_bc_matrix");

Maize.xu48FM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu48_raw/raw_feature_bc_matrix");


teosenti.xu53FM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu53_raw/raw_feature_bc_matrix");

teosenti.xu55FM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu55_raw/raw_feature_bc_matrix");


## Back to integration folder 

#...Desktop/Seurat/Maize-vs-teosinte. 

# My example below
setwd("C:/Lab_Xu/Maize-vs-teosinte")

## Convert each feature-barcode matrix to a Seurat object.

Maize.xu43FM <- CreateSeuratObject(counts = Maize.xu43FM.data, min.cells=3, min.features=200, project="Maize.xu43FM");


Maize.xu48FM <- CreateSeuratObject(counts = Maize.xu48FM.data, min.cells=3, min.features=200, project="Maize.xu48FM");


teosenti.xu53FM <- CreateSeuratObject(counts = teosenti.xu53FM.data, min.cells=3, min.features=200, project="teosenti.xu53FM");


teosenti.xu55FM <- CreateSeuratObject(counts = teosenti.xu55FM.data, min.cells=3, min.features=200, project="teosenti.xu55FM");



## View metadata to see the data information. 
## origin.ident
## nCount_RNA
## nFeature_RNA

# View(Maize.xu43FM@meta.data)
# View(Maize.xu48FM@meta.data)
# 
# View(teosenti.xu53FM@meta.data)
# View(teosenti.xu55FM@meta.data)

#View(teosinte.4@meta.data)

#2
# VlnPlot(Maize.xu43FM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
# VlnPlot(Maize.xu48FM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
# VlnPlot(teosenti.xu53FM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
# VlnPlot(teosenti.xu55FM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")

Maize.xu43FM <- subset(Maize.xu43FM, subset = nCount_RNA > 5000 & nCount_RNA < 150000 & nFeature_RNA>1000 &nFeature_RNA<15000)

Maize.xu48FM <- subset(Maize.xu48FM, subset = nCount_RNA > 5000 & nCount_RNA < 150000 & nFeature_RNA>1000 &nFeature_RNA<15000)



teosenti.xu53FM <- subset(teosenti.xu53FM, subset = nCount_RNA > 5000 & nCount_RNA < 150000 & nFeature_RNA>1000 &nFeature_RNA<15000)

teosenti.xu55FM <- subset(teosenti.xu55FM, subset = nCount_RNA > 5000 & nCount_RNA < 150000 & nFeature_RNA>1000 &nFeature_RNA<15000)
# VlnPlot(Maize.xu43FM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
# VlnPlot(Maize.xu48FM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
# VlnPlot(teosenti.xu53FM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
# VlnPlot(teosenti.xu55FM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")

Maize.xu43FM@meta.data$source <- "Maize"
Maize.xu48FM@meta.data$source <- "Maize"


teosenti.xu53FM@meta.data$source <- "teosinte"
teosenti.xu55FM@meta.data$source <- "teosinte"

## View metadata to see the data information. 
## origin.ident
## nCount_RNA
## nFeature_RNA
## source

# View(Maize.xu43FM@meta.data)
# View(Maize.xu48FM@meta.data)
# 
# View(teosenti.xu53FM@meta.data)
# View(teosenti.xu55FM@meta.data)
#View(teosinte.4@meta.data)




# 4. merge all samples into one object (obj) for integrated analysis



Maize_vs_teosinte.merge <- merge(x= Maize.xu43FM, y=list(Maize.xu48FM, teosenti.xu53FM, teosenti.xu55FM))
View(Maize_vs_teosinte.merge@meta.data)



# 5. Perform analysis without integration

# We can first analyze the dataset without integration. The resulting clusters are defined both by cell type and treatment condition, which may create challenges for downstream analysis.

# Option-1 
# run standard LongNormalization anlaysis workflow
Maize_vs_teosinte <- NormalizeData(Maize_vs_teosinte.merge)
Maize_vs_teosinte <- FindVariableFeatures(Maize_vs_teosinte)
Maize_vs_teosinte <- ScaleData(Maize_vs_teosinte)


# Option-2 
# we can also perform integration using sctransform-normalized data as we did for individual dataset.

# options(future.globals.maxSize = 3e+09)
# Maize_vs_teosinte <- SCTransform(Maize_vs_teosinte)

###




# 6. Perform linear dimensional reduction



# RunPCA

Maize_vs_teosinte <- RunPCA(Maize_vs_teosinte)

# Determine the ‘dimensionality’ of the dataset


## We will use a quantitative approach based on Elbow plot: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function). 
## https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html

## We first display the firt 50 PCs on the ElbowPlot. 

ElbowPlot(Maize_vs_teosinte, ndims = 50)

## Based on ElbowPlot(object, ndims = 50) plot, we could roughly determine the majority of the variation by where the elbow occurs (touches the ground). While this gives us a good rough idea of the number of PCs needed to be included, a more quantitative approach may be a bit more reliable. We can calculate where the principal components start to elbow by taking the larger (should be smaller?) value of:

## (1) The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
## (2) The point where the percent change in variation between the consecutive PCs is less than 0.1%.

## We will start by calculating the first metriE:
# Determine percent of variation associated with each PC
pct <- Maize_vs_teosinte[["pca"]]@stdev / sum(Maize_vs_teosinte[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

## [1] result will be displayed (i.e. 43)

## The first metric returns PC (i.e. 43) as the PC matching these requirements. Let’s check the second metric, which identifies the PC where the percent change in variation between consecutive PCs is less than 0.1%:

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

## [1]  result will be displayed (i.e. 17)

## This second metric returns PC (i.e. 17). Usually, we would choose the minimum of these two metrics as the PCs covering the majority of the variation in the data.

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs




# 7. Cluster the cells

Maize_vs_teosinte <- IntegrateLayers(
    object = Maize_vs_teosinte, method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE
)

# Now let's try one by one, and then compare the cluster results.

# CCAIntegration
dims_used_list <- list(1:17, 1:25, 1:30)
resolution_used_list <- c(1, 1.25, 1.5)
combined_plots_list <- list()
i <- 1  
marker_genes <- c("GRMZM2G017087", "GRMZM2G001289", "GRMZM2G083725", "GRMZM2G372364", "GRMZM2G165836", "GRMZM2G047448", "GRMZM2G466532", "GRMZM2G354253", "GRMZM2G073671", "GRMZM2G310115", "GRMZM2G005619", "GRMZM2G305046", "GRMZM2G071959", "GRMZM2G119071", "GRMZM2G102161", "GRMZM2G159397", "GRMZM2G307119", "GRMZM2G097989", "GRMZM2G014729", "GRMZM2G003927", "GRMZM2G397518", "GRMZM2G072274", "GRMZM2G529859", "GRMZM2G005353", "GRMZM2G088309", "GRMZM2G102218", "GRMZM2G159399", "GRMZM2G043854", "GRMZM2G072342", "GRMZM2G526668", "GRMZM2G176141", "GRMZM2G052544", "GRMZM2G140694", "GRMZM2G589696", "GRMZM2G017470", "GRMZM2G144172", "GRMZM2G132794", "GRMZM2G345700", "GRMZM2G361652", "GRMZM2G430849", "GRMZM2G039074")
for (dims_used in dims_used_list){
  dims_range_pdf <- paste0(min(dims_used), "_", max(dims_used))
  for (resolution_used in resolution_used_list){
    
    Maize_vs_teosinte <- FindNeighbors(Maize_vs_teosinte, dims = dims_used, reduction = "integrated.cca")
    Maize_vs_teosinte <- FindClusters(Maize_vs_teosinte, resolution = resolution_used, cluster.name = "cca_clusters")
    Maize_vs_teosinte <- RunUMAP(Maize_vs_teosinte, dims = dims_used, reduction = "integrated.cca", reduction.name = "umap.cca")
    title_text <- paste0("FM - Dims: ", dims_range_pdf, ", Res: ", resolution_used)
    
    p4 <- DimPlot(
      Maize_vs_teosinte,
      reduction = "umap.cca",
      label = TRUE,
      pt.size = 0.5,
      label.size = 8
    ) + ggtitle(title_text)
    
    p5 <- DimPlot(
      Maize_vs_teosinte,
      reduction = "umap.cca",
      group.by = "orig.ident",
      pt.size = 0.5
    ) + ggtitle(title_text)
    
    p6 <- DimPlot(
      Maize_vs_teosinte,
      reduction = "umap.cca",
      group.by = "source",
      pt.size = 0.5
    ) + ggtitle(title_text)
    
    # Generate the violin plot
    vln <- VlnPlot(
      object = Maize_vs_teosinte,
      features = c("nFeature_RNA", "nCount_RNA"),
      ncol = 2
    ) + ggtitle(title_text)
    
    # Arrange the four plots into a 2x2 grid
    combined_plot <- grid.arrange(p4, p5, p6, vln, ncol = 4)
    # Alternatively, using cowplot:
    # combined_plot <- plot_grid(p4, p5, p6, vln, ncol = 2)
    # Or using patchwork:
    # combined_plot <- (p4 + p5) / (p6 + vln)
    
    # Store the combined plot
    combined_plots_list[[i]] <- combined_plot
    names(combined_plots_list)[i] <- paste0("Dims_", dims_used, "_Res_", resolution_used)
    i <- i + 1 
  }
  markercca<-FeaturePlot(Maize_vs_teosinte, features = marker_genes, cols = c("light blue", "red"), reduction = "umap.cca", pt.size = 1.5)
  pdf_filename <- paste0("FM_markergene_Plots_Dims_", dims_range_pdf, ".pdf")
  pdf(pdf_filename, width = 40, height = 135)
  grid.newpage()
  grid.draw(markercca)
  dev.off()
}
final_plot <- arrangeGrob(grobs = combined_plots_list, ncol = 1)

# Open a PDF device with the desired dimensions
pdf("FM_All_Combined_Plots.pdf", width = 45.6, height = 135)

# Draw the final combined plot
grid.newpage()
grid.draw(final_plot)

# Close the PDF device
dev.off()

## Further discussion for Option-2 scTransformation. 
# Option-2 
# If we perform integration using sctransform-normalized data as we did for individual dataset.we need to modify the code to replace longNormalization as example below. 

#options(future.globals.maxSize = 3e+09)
#Maize_vs_teosinte <- SCTransform(Maize_vs_teosinte)
#Maize_vs_teosinte <- RunPCA(Maize_vs_teosinte, npcs = 30, verbose = F)


#Maize_vs_teosinte <- IntegrateLayers(
#  object = obj,
#  method = RPCAIntegration,
#  normalization.method = "SCT",
#  verbose = F
#)
#Maize_vs_teosinte <- FindNeighbors(Maize_vs_teosinte, dims = 1:17, reduction = "integrated.dr")
#Maize_vs_teosinte <- FindClusters(Maize_vs_teosinte, resolution = 2)


## Check gene expression. 

# Define marker genes
marker_genes <- c("GRMZM2G017087", "GRMZM2G001289", "GRMZM2G083725", "GRMZM2G372364", "GRMZM2G165836", "GRMZM2G047448", "GRMZM2G466532", "GRMZM2G354253", "GRMZM2G073671", "GRMZM2G310115", "GRMZM2G005619", "GRMZM2G305046", "GRMZM2G071959", "GRMZM2G119071", "GRMZM2G102161", "GRMZM2G159397", "GRMZM2G307119", "GRMZM2G097989", "GRMZM2G014729", "GRMZM2G003927", "GRMZM2G397518", "GRMZM2G072274", "GRMZM2G529859", "GRMZM2G005353", "GRMZM2G088309", "GRMZM2G102218", "GRMZM2G159399", "GRMZM2G043854", "GRMZM2G072342", "GRMZM2G526668", "GRMZM2G176141", "GRMZM2G052544", "GRMZM2G140694", "GRMZM2G589696", "GRMZM2G017470", "GRMZM2G144172", "GRMZM2G132794", "GRMZM2G345700", "GRMZM2G361652", "GRMZM2G430849", "GRMZM2G039074")

# umap.cca
markercca<-FeaturePlot(Maize_vs_teosinte, features = marker_genes, cols = c("light blue", "red"), reduction = "umap.cca")
pdf("FM_markergene_Plots.pdf", width = 40, height = 135)
markercca
# Draw the final combined plot
grid.newpage()
grid.draw(markercca)

# Close the PDF device
dev.off()

## Compare individual gene for all integration method side by side

# Define marker genes
marker_single_genes <- c("GST4")



marker_single.cca <- FeaturePlot(Maize_vs_teosinte, features = marker_single_genes, cols = c("light blue", "red"), reduction = "umap.cca")
marker_single.cca

###############################################################################################################################################################################
# 11. Finding differentially expressed features/genes between Maize and teosinte treatment in different cell types. 




# Once integrative analysis is complete, you can rejoin the layers - which collapses the individual datasets together and recreates the original counts and data layers. You will need to do this before performing any differential expression analysis. However, you can always resplit the layers in case you would like to reperform integrative analysis.

# assign a new object name so it won't change the original one. 

Maize_vs_teosinte_deg <- Maize_vs_teosinte



## Without Joinlayers, we still have individual dataset in RNA/layers, run the code below: 

Layers(Maize_vs_teosinte_deg[["RNA"]])

# [1] "counts.Maize.xu45SM" "counts.Maize.xu47SM" "counts.Maize.3" "counts.Maize.4" "counts.teosenti.xu52SM"   
#  [6] "counts.teosenti.xu56SM"    "counts.teosenti.xu58SM"    "counts.teosinte.4"    "data.Maize.xu45SM"   "data.Maize.xu47SM"  
# [11] "data.Maize.3"   "data.Maize.4"   "data.teosenti.xu52SM"      "data.teosenti.xu56SM"      "data.teosenti.xu58SM"     
# [16] "data.teosinte.4"      "scale.data"    




# rejoin the layers

Maize_vs_teosinte_deg <- JoinLayers(Maize_vs_teosinte_deg)



# after Joinlayers, we only have in RNA/layers, run the code below: 

Layers(Maize_vs_teosinte_deg[["RNA"]])

#[1] "data"       "counts"     "scale.data"



## Find DEG for certain cluster between the two conditions (Maize vs teosinte treatment). 

# Specify the reduction method
reduction_method <- "umap.cca"

# Set the identity class to use 'cca_clusters' and assign to new object
Maize_vs_teosinte_deg <- SetIdent(Maize_vs_teosinte_deg, value = "cca_clusters")

# Combine cluster identities with treatment information
Maize_vs_teosinte_deg$cluster_treatment <- paste(Idents(Maize_vs_teosinte_deg), Maize_vs_teosinte_deg$source, sep = "_")
for (cluster_num in 0:26) {
  
  Idents(Maize_vs_teosinte_deg) <- "cluster_treatment"
  ident1 <- paste0(cluster_num, "_Maize")     
  ident2 <- paste0(cluster_num, "_teosinte") 
  degs <- FindMarkers(
    object = Maize_vs_teosinte_deg,
    ident.1 = ident1,
    ident.2 = ident2,
    assay = "RNA",
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )
  head(degs, n = 15)
  
  # Generate the filename dynamically based on the reduction method
  filename <- paste0("DEGs_FM_Cluster", cluster_num, "_Maize_vs_teosinte_", reduction_method, ".xls")
  
  # Export the DEGs to a tab-separated text file
  write.table(degs, file = filename, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

}
saveRDS(Maize_vs_teosinte, file = "Maize_vs_teosinte_LogNormalizeFM.rds")
saveRDS(Maize_vs_teosinte_deg, file = "Maize_vs_teosinte_LogNormalizeFM_Deg.rds")

# Or, you may save more multiple object in one rds file. 
# Combine the objects into a list
Maize_vs_teosinte_all <- list(
  LogNormalizeFM = Maize_vs_teosinte,
  LogNormalizeFM_Filtered = Maize_vs_teosinte_filtered,
  LogNormalizeFM_Deg = Maize_vs_teosinte_deg
)

# Save the list to a single RDS file
saveRDS(Maize_vs_teosinte_all, file = "Maize_vs_teosinte_FM_all.rds")

#Negative is Teosenti
