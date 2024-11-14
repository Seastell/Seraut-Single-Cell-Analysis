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
setwd("C:/Lab_Xu/Data/xu29_raw")

setwd("C:/Lab_Xu/Data/xu32_raw")

setwd("C:/Lab_Xu/Data/xu49_raw")

setwd("C:/Lab_Xu/Data/xu51_raw")

Maize.xu29IM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu29_raw/raw_feature_bc_matrix");

Maize.xu32IM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu32_raw/raw_feature_bc_matrix");

#Maize.xu44IM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu44_raw/raw_feature_bc_matrix");

#Maize.xu46IM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu46_raw/raw_feature_bc_matrix");

teosenti.xu49IM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu49_raw/raw_feature_bc_matrix");

teosenti.xu51IM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu51_raw/raw_feature_bc_matrix");

#teosenti.xu57IM.data <- Read10X(data.dir = "C:/Lab_Xu/Data/xu57_raw/raw_feature_bc_matrix");


## Back to integration folder 

#...Desktop/Seurat/Maize-vs-teosinte. 

# My example below
setwd("C:/Lab_Xu/Maize-vs-teosinte")

## Convert each feature-barcode matrix to a Seurat object.

Maize.xu29IM <- CreateSeuratObject(counts = Maize.xu29IM.data, min.cells=3, min.features=200, project="Maize.xu29IM");


Maize.xu32IM <- CreateSeuratObject(counts = Maize.xu32IM.data, min.cells=3, min.features=200, project="Maize.xu32IM");


#Maize.xu44IM <- CreateSeuratObject(counts = Maize.xu44IM.data, min.cells=3, min.features=200, project="Maize.xu44IM");


#Maize.xu46IM <- CreateSeuratObject(counts = Maize.xu46IM.data, min.cells=3, min.features=200, project="Maize.xu46IM");


teosenti.xu49IM <- CreateSeuratObject(counts = teosenti.xu49IM.data, min.cells=3, min.features=200, project="teosenti.xu49IM");


teosenti.xu51IM <- CreateSeuratObject(counts = teosenti.xu51IM.data, min.cells=3, min.features=200, project="teosenti.xu51IM");


#teosenti.xu57IM <- CreateSeuratObject(counts = teosenti.xu57IM.data, min.cells=3, min.features=200, project="teosenti.xu57IM");

## View metadata to see the data information. 
## origin.ident
## nCount_RNA
## nFeature_RNA

#View(Maize.xu29IM@meta.data)
#View(Maize.xu32IM@meta.data)
#View(Maize.xu44IM@meta.data)
#View(Maize.xu46IM@meta.data)
#View(teosenti.xu49IM@meta.data)
#View(teosenti.xu51IM@meta.data)
#View(teosenti.xu57IM@meta.data)
#View(teosinte.4@meta.data)

#2
VlnPlot(Maize.xu29IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
VlnPlot(Maize.xu32IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
#VlnPlot(Maize.xu44IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
#VlnPlot(Maize.xu46IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
VlnPlot(teosenti.xu49IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
VlnPlot(teosenti.xu51IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
#VlnPlot(teosenti.xu57IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "Before Filtering")
#VlnPlot(Maize.xu29IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)
Maize.xu29IM <- subset(Maize.xu29IM, subset = nCount_RNA > 5000 & nCount_RNA < 110000 & nFeature_RNA>1000 & nFeature_RNA<13000)

Maize.xu32IM <- subset(Maize.xu32IM, subset = nCount_RNA > 5000 & nCount_RNA < 100000 & nFeature_RNA>1000 & nFeature_RNA<13000)

#Maize.xu44IM <- subset(Maize.xu44IM, subset = nCount_RNA > 10000 & nCount_RNA < 100000 & nFeature_RNA>1000 & nFeature_RNA<15000)

#Maize.xu46IM <- subset(Maize.xu46IM, subset = nCount_RNA > 10000 & nCount_RNA < 150000 & nFeature_RNA>500 & nFeature_RNA<15000)

teosenti.xu49IM <- subset(teosenti.xu49IM, subset = nCount_RNA > 5000 & nCount_RNA < 40000 & nFeature_RNA>500 & nFeature_RNA<10000)

teosenti.xu51IM <- subset(teosenti.xu51IM, subset = nCount_RNA > 500 & nCount_RNA < 15000 & nFeature_RNA>100 & nFeature_RNA<6000)

#teosenti.xu57IM <- subset(teosenti.xu57IM, subset = nCount_RNA > 50 & nCount_RNA < 800 & nFeature_RNA>50 & nFeature_RNA<650)
#VlnPlot(Maize.xu29IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
#VlnPlot(Maize.xu32IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
#VlnPlot(Maize.xu44IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
#VlnPlot(Maize.xu46IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
#VlnPlot(teosenti.xu49IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
#VlnPlot(teosenti.xu51IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
#VlnPlot(teosenti.xu57IM, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)+plot_annotation(title = "After Filtering")
Maize.xu29IM@meta.data$source <- "Maize"
Maize.xu32IM@meta.data$source <- "Maize"
#Maize.xu44IM@meta.data$source <- "Maize"
#Maize.xu46IM@meta.data$source <- "Maize"

teosenti.xu49IM@meta.data$source <- "teosinte"
teosenti.xu51IM@meta.data$source <- "teosinte"
#teosenti.xu57IM@meta.data$source <- "teosinte"
#teosinte.4@meta.data$source <- "teosinte"
## View metadata to see the data information. 
## origin.ident
## nCount_RNA
## nFeature_RNA
## source

#View(Maize.xu29IM@meta.data)
#View(Maize.xu32IM@meta.data)
#View(Maize.xu44IM@meta.data)
#View(Maize.xu46IM@meta.data)
#View(teosenti.xu49IM@meta.data)
#View(teosenti.xu51IM@meta.data)
#View(teosenti.xu57IM@meta.data)
#View(teosinte.4@meta.data)

# 4. merge all samples into one object (obj) for integrated analysis

Maize_vs_teosinte.merge <- merge(x= Maize.xu29IM, y=list(Maize.xu32IM, teosenti.xu49IM, teosenti.xu51IM))
#View(Maize_vs_teosinte.merge@meta.data)



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


# Cluster the cell. 

## View metadata to see the data information. 
## origin.ident
## nCount_RNA
## nFeature_RNA
## source
## unintegrated_clusters
## seurat_clusters

#View (Maize_vs_teosinte@meta.data)

# CCAIntegration
###########################################################################################################

Maize_vs_teosinte <- IntegrateLayers(
  object = Maize_vs_teosinte, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

dims_used_list <- list(1:16, 1:25, 1:30)
resolution_used_list <- c(1, 1.25, 1.5)
combined_plots_list <- list()
i <- 1  
marker_genes <- c("GRMZM2G017087", "GRMZM2G001289", "GRMZM2G083725", "GRMZM2G372364", "GRMZM2G165836", "GRMZM2G047448", "GRMZM2G466532", "GRMZM2G354253", "GRMZM2G073671", "GRMZM2G310115", "GRMZM2G005619", "GRMZM2G305046", "GRMZM2G071959", "GRMZM2G119071", "GRMZM2G102161", "GRMZM2G159397", "GRMZM2G307119", "GRMZM2G097989", "GRMZM2G014729", "GRMZM2G003927", "GRMZM2G397518", "GRMZM2G072274", "GRMZM2G529859", "GRMZM2G005353", "GRMZM2G088309", "GRMZM2G102218", "GRMZM2G159399", "GRMZM2G043854", "GRMZM2G072342", "GRMZM2G526668", "GRMZM2G176141", "GRMZM2G052544", "GRMZM2G140694", "GRMZM2G589696", "GRMZM2G017470", "GRMZM2G144172", "GRMZM2G132794", "GRMZM2G345700", "GRMZM2G361652", "GRMZM2G430849", "GRMZM2G039074")
for (dims_used in dims_used_list){
  dims_range_pdf <- paste0(min(dims_used), "_", max(dims_used))
  for (resolution_used in resolution_used_list){
    
    Maize_vs_teosinte <- FindNeighbors(Maize_vs_teosinte, dims = dims_used, reduction = "integrated.cca")
    Maize_vs_teosinte <- FindClusters(Maize_vs_teosinte, resolution = resolution_used, cluster.name = "cca_clusters")
    # 8. Run non-linear dimensional reduction (UMAP/tSNE)
    Maize_vs_teosinte <- RunUMAP(Maize_vs_teosinte, dims = dims_used, reduction = "integrated.cca", reduction.name = "umap.cca")
    title_text <- paste0("IM - Dims: ", dims_range_pdf, ", Res: ", resolution_used)
    
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
  pdf_filename <- paste0("IM_markergene_Plots_Dims_", dims_range_pdf, ".pdf")
  pdf(pdf_filename, width = 40, height = 135)
  grid.newpage()
  grid.draw(markercca)
  dev.off()
}

final_plot <- arrangeGrob(grobs = combined_plots_list, ncol = 1)

pdf("IM_All_Combined_Plots.pdf", width = 45.6, height = 135)

# Draw the final combined plot
grid.newpage()
grid.draw(final_plot)

# Close the PDF device
dev.off()


##############################################################################################################################

Maize_vs_teosinte_filtered <- subset(Maize_vs_teosinte, subset = cca_clusters != 0)
Maize_vs_teosinte_filtered <- IntegrateLayers(
  object = Maize_vs_teosinte_filtered, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
Maize_vs_teosinte_filtered <- FindNeighbors(Maize_vs_teosinte_filtered, dims = 1:30, reduction = "integrated.cca")
Maize_vs_teosinte_filtered <- FindClusters(Maize_vs_teosinte_filtered, resolution = 1.5, cluster.name = "cca_clusters")

Maize_vs_teosinte_filtered <- RunUMAP(Maize_vs_teosinte_filtered, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap.cca")

p4 <- DimPlot(Maize_vs_teosinte_filtered, reduction = "umap.cca", label = TRUE, pt.size = 1, label.size = 8)

p5 <- DimPlot(Maize_vs_teosinte_filtered, reduction = "umap.cca", group.by = "orig.ident", pt.size = 1)

p6 <- DimPlot(Maize_vs_teosinte_filtered, reduction = "umap.cca", group.by = "source", pt.size = 1)
vln <- VlnPlot(
  object = Maize_vs_teosinte_filtered,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2
)
combined_plot_filtered <- grid.arrange(p4, p5, p6, vln, ncol = 4)
pdf("IM_All_Combined_Filtered_Plots.pdf", width = 45.6, height = 15)
grid.newpage()
grid.draw(combined_plot_filtered)
dev.off()

## Check gene expression. 

# Define marker genes
marker_genes <- c("GRMZM2G017087", "GRMZM2G001289", "GRMZM2G083725", "GRMZM2G372364", "GRMZM2G165836", "GRMZM2G047448", "GRMZM2G466532", "GRMZM2G354253", "GRMZM2G073671", "GRMZM2G310115", "GRMZM2G005619", "GRMZM2G305046", "GRMZM2G071959", "GRMZM2G119071", "GRMZM2G102161", "GRMZM2G159397", "GRMZM2G307119", "GRMZM2G097989", "GRMZM2G014729", "GRMZM2G003927", "GRMZM2G397518", "GRMZM2G072274", "GRMZM2G529859", "GRMZM2G005353", "GRMZM2G088309", "GRMZM2G102218", "GRMZM2G159399", "GRMZM2G043854", "GRMZM2G072342", "GRMZM2G526668", "GRMZM2G176141", "GRMZM2G052544", "GRMZM2G140694", "GRMZM2G589696", "GRMZM2G017470", "GRMZM2G144172", "GRMZM2G132794", "GRMZM2G345700", "GRMZM2G361652", "GRMZM2G430849", "GRMZM2G039074")

# umap.cca_filtered
markercca_filtered<-FeaturePlot(Maize_vs_teosinte_filtered, features = marker_genes, cols = c("light blue", "red"), reduction = "umap.cca",pt.size=1.5)
pdf("IM_markergene_Plots_filtered.pdf", width = 40, height = 135)
# Draw the final combined plot
grid.newpage()
grid.draw(markercca_filtered)
dev.off()

###############################################################################################################################################################################
# 11. Finding differentially expressed features/genes between Maize and teosinte treatment in different cell types. 




# Once integrative analysis is complete, you can rejoin the layers - which collapses the individual datasets together and recreates the original counts and data layers. You will need to do this before performing any differential expression analysis. However, you can always resplit the layers in case you would like to reperform integrative analysis.

# assign a new object name so it won't change the original one. 

Maize_vs_teosinte_deg <- Maize_vs_teosinte_filtered



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
for (cluster_num in 0:25) {
  
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
  filename <- paste0("DEGs_IM_Cluster", cluster_num, "_Maize_vs_teosinte_", reduction_method, ".xls")
  
  # Export the DEGs to a tab-separated text file
  write.table(degs, file = filename, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

}
saveRDS(Maize_vs_teosinte, file = "Maize_vs_teosinte_LogNormalizeIM.rds")
saveRDS(Maize_vs_teosinte_filtered, file = "Maize_vs_teosinte_LogNormalizeIM_Filtered.rds")
saveRDS(Maize_vs_teosinte_deg, file = "Maize_vs_teosinte_LogNormalizeIM_Deg.rds")

# Or, you may save more multiple object in one rds file. 
# Combine the objects into a list
Maize_vs_teosinte_all <- list(
  LogNormalizeIM = Maize_vs_teosinte,
  LogNormalizeIM_Filtered = Maize_vs_teosinte_filtered,
  LogNormalizeIM_Deg = Maize_vs_teosinte_deg
)

# Save the list to a single RDS file
saveRDS(Maize_vs_teosinte_all, file = "Maize_vs_teosinte_IM_all.rds")


