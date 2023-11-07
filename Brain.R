# Richard Mustaklem - Spring 2023
# Student of Dr. Kim Jungsu and apprentice of Dom Acri
# scRNA-seq Seurat Tutorial: The Basics


#1 - Load libraries
  library(dplyr)
  library(Seurat)
  library(patchwork)
 
#2 - Download the data
  # These instructions are out of date
  # Go to https://www.10xgenomics.com/
  # Click on Resources -> Datasets
  # Click on Single Cell Gene Expression (left side)
  # Click on Nuclei Isolation for Single Cell Expression
  # Click on 5k Adult Mouse Brain Nuclei Isolated with Chromium Nuclei Isolation Kit
  # Click on Feature / cell matrix HDF5 (filtered)
  # Extract the "filtered_feature_bc_matrix"

#3 - Load dataset
  #Check if the file is in your working directory
  dir()
  #Change your working directory to location of files
  setwd("C:/Users/Richie/Documents/")
  #Load in 10x Genomics data
  brain.data <- Read10X("filtered_feature_bc_matrix")
  #Create RNA object, where the counts
  brain <- CreateSeuratObject(counts = brain.data, project = "MBD", min.cells = 3, min.features = 200)
  brain # Features, samples, 1 assay
  head(brain) # We can look at different column data
  
#4 -  QC and selecting cells for further analysis
  # Each species has their own identification/patterns or 
  # rules to help guide researchers on what they work on.
  # Capitalized variables, patterns, objects, etc, are representative of HUMAN data. 
  # For example, this data is mouse data, so this data uses lower case "^mt" 
  # rather than "^MT" which was used in the Seurat data set (human)
  brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-") # Add percent.mt column, showing genes with ^mt in the name
  head(brain) # We can look at different column data, notice theres now percent.mt column
  VlnPlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Display these column types 
  
  # The mitochondria value is still considered and this can lead to a lot of noise.
  # Percent.mt can be tricky to cut off, because a lot of times there are
  # dead cells that still get captured in data.
  # Around 15-20 percent can be a decent spot to cut off. Each scenario is case by case, 
  # however some use the 5-10% range to ensure live cells. 
  # There can be lots of money lost with just throwing out large amount of data.
  brain <- subset(brain, subset = nFeature_RNA < 6000 & nCount_RNA < 30000 & percent.mt < 15)
  VlnPlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#5 -  Normalize the Data
  # Highly variable gene shown pre-normalization
  VlnPlot(brain,features="Hbb-bs") # Try Inpp5d as well, a marker for microglia
  # Normalization
  brain <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000) # Scale factor is the magnitude of your normalization. 10000 is recommend by Satija
  # Highly variable gene shown pre-normalization
  VlnPlot(brain,features="Hbb-bs") # Look at the difference after normalization
  
#6 - Identification of highly variable features
  # Select only the top 2000 genes 
  brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000) # Identifies features that are outliers on a 'mean variability plot'.
  
  # The code below is just for fun, not really required for analysis.
  top10 <- head(VariableFeatures(brain), 10) # Find the top ten variable genes, display them
  top10
  plot1 <- VariableFeaturePlot(brain)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
  plot1 + plot2
  plot1
  plot2
  
#7 - Scaling the data
  all.genes <- rownames(brain)
  brain <- ScaleData(brain, features = all.genes)
  
#8 -Perform linear dimensional reduction
  brain <- RunPCA(brain, features = VariableFeatures(object = brain))
  VizDimLoadings(brain, dims = 1:2, reduction = "pca")
  DimPlot(brain, reduction = "pca")
  
  # There is clear association though in the first 7 or 8 PC
  # Question: how do we make sure that everything is going correctly
  DimHeatmap(brain, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(brain, dims = 1:15, cells = 500, balanced = TRUE)
  
#9 - Determine the dimensionality of the dataset
  brain <- JackStraw(brain, num.replicate = 100)
  brain <- ScoreJackStraw(brain, dims = 1:20)
  ElbowPlot(brain,ndims = 50)

#10 - # Cluster the cells
  # Use elbow plot to determine how many dimensions to look at, 
  # Here we can go to 1:25 but even 1:40 works
  brain <- FindNeighbors(brain, dims = 1:25)
  brain <- FindClusters(brain, resolution = 0.5) # Think of resolutions like microscopes. You can use different resolutions
                                                 # to change how "zoomed in" you want to be

  brain@active.ident # Looking at the current amount of clusters we are actually using
  Idents(brain) <- "RNA_snn_res.0.6" # Can switch between any of the resolutions with this command
  
#11 - Run non-linear dimensional reduction (UMAP/tSNE)
  brain <- RunUMAP(brain, dims = 1:25)
  DimPlot(brain, reduction = "umap", label = TRUE)
  # Example below on looking at one specific gene, Hbb-bs
  FeaturePlot(brain, features = c("Hbb-bs")) # Let's look at Inpp5d
  
  # Save your R file to save computation time in the future
  # saveRDS(brain,file = "file_name.rds")

#12 - # Finding differentially expressed features
  # Get all genes that are positive + assign markers
  neuron.markers <- FindAllMarkers(brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  head(neuron.markers)
  # How to see the 1st, second, third cluster?
  head(subset(neuron.markers, cluster =="11"))
  # Then go to https://panglaodb.se/samples.html and manually look up genes, 
  # Look at each cluster, significant genes (p-value), 
  # make 25 different names (0:24)
  # Rename clusters ID
  new.cluster.ids <- c("Oligodendrocyte", "Neurons", "Unknown", "Astrocytes", "Neurons", " Neurons",
                       "Germ", "Macrophages", "Germ", "Enterocytes", "Neurons",
                       "Osteoblast", "OPC", "Unknown", "Fibroblast", "Unknown", "Fibroblast",
                       "Endothelial", "Fibroblast", "Fibroblast", "Neurons", "Unknown", "Endothelial",
                       "T Cell")
  names(new.cluster.ids) <- levels(brain)
  brain <- RenameIdents(brain, new.cluster.ids)
  # Plot the Umap but with the annotated clusters
  DimPlot(brain, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  
  