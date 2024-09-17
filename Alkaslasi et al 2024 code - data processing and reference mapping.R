
---
title: "Alkaslasi et al. 2024"
output: html_document
date: "2024-07-20"
---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

setwd("/Users/alkaslasimr/Documents/Sequencing/TBISeq/")
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
#devtools::install_github(repo = "samuel-marsh/scCustomize")
library(scCustomize)
library(clipr)

```


```{r load_data}
## Load un-processed data
load("~/Documents/Sequencing/TBISeq/MA02_021__Atf3CreER_M7dpi/Atf3CreER_7dpi_M.data.robj") # Atf3-CreER males
load("~/Documents/Sequencing/TBISeq/MA02_004 Atf3-CreER/Atf3CreER_7dpi_F.data.robj") # Atf3-CreER females

atf3_m <-CreateSeuratObject(Atf3CreER_7dpi_M.data)
atf3_f <-CreateSeuratObject(Atf3CreER_7dpi_F.data)
```

```{r process_data_QC, echo=FALSE}
## Load un-processed data
obj.list <- list(atf3_m = atf3_m,
                 atf3_f = atf3_f)

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- PercentageFeatureSet(object = x, pattern = "^mt-", col.name = "percent.mt")
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- NormalizeData(x)
  x <- SCTransform(object = x, vars.to.regress = c("percent.mt"), verbose = FALSE)
  #x <- subset(x, subset = percent.mt < 5)
  x <- RunPCA(object = x, npcs = 50)
  x <- RunUMAP(object = x, dims = 1:30, verbose = FALSE)
  x <- FindNeighbors(object = x, dims = 1:30, verbose = FALSE)
  x <- FindClusters(object = x, verbose = FALSE)
})

atf3_m_nofilt <- obj.list[[1]]
atf3_f_nofilt <- obj.list[[2]]

VlnPlot(atf3_m_nofilt, group.by = 'orig.ident', features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(atf3_f_nofilt, group.by = 'orig.ident', features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(atf3_m_nofilt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(atf3_m_nofilt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(atf3_f_nofilt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(atf3_f_nofilt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
```

```{r process_data, echo=FALSE}

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- PercentageFeatureSet(object = x, pattern = "^mt-", col.name = "percent.mt")
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- NormalizeData(x)
  x <- SCTransform(object = x, vars.to.regress = c("percent.mt"), verbose = FALSE)
  x <- subset(x, subset = percent.mt < 5)
  x <- RunPCA(object = x, npcs = 50)
  x <- RunUMAP(object = x, dims = 1:30, verbose = FALSE)
  x <- FindNeighbors(object = x, dims = 1:30, verbose = FALSE)
  x <- FindClusters(object = x, verbose = FALSE)
})

atf3_m <- obj.list[[1]]
atf3_f <- obj.list[[2]]


atf3_m$group <- "Atf3"
atf3_m$expt <- "Atf3-CreER_M"
atf3_m$sex <- "M"

atf3_f$group <- "Atf3"
atf3_f$expt <- "Atf3-CreER_F"
atf3_f$sex <- "F"

DimPlot(atf3_m) + ggtitle("Atf3-CreER Males")
DimPlot(atf3_f) + ggtitle("Atf3-CreER Females")

```


```{r integrate_atf3, echo=FALSE}
## Integrate Atf3 samples

obj.list2 <- c(atf3_m = atf3_m,
               atf3_f = atf3_f)

features <- SelectIntegrationFeatures(object.list = obj.list2)
obj.list2 <- PrepSCTIntegration(object.list = obj.list2, anchor.features = features)

control.anchors <- FindIntegrationAnchors(object.list = obj.list2, normalization.method = "SCT",
                                          anchor.features = features)
int_obj <- IntegrateData(anchorset = control.anchors, normalization.method = "SCT")
int_obj <- RunPCA(int_obj, verbose = FALSE)
int_obj <- RunUMAP(int_obj, reduction = "pca", dims = 1:30)
int_obj <- FindNeighbors(object = int_obj, reduction = "pca", dims = 1:30)
int_obj <- FindClusters(int_obj, resolution = 0.8)

atf3creer <- int_obj

DimPlot(atf3creer, split.by = 'expt')

```


```{r process_refdata, echo=FALSE}

refdata <- readRDS("~/Documents/Sequencing/TBISeq/addressing reviewer comments/allen_mop_2020.rds") # data downloaded directly from https://seurat.nygenome.org/azimuth/demo_datasets/allen_mop_2020.rds
#newref <- readRDS(file.choose()) # data downloaded directly from https://zenodo.org/records/4546935

# Try just mapping the subset data from azimuth to the azimuth ref to make sure it has those predicted subclass
projected.umap <- readRDS("~/Documents/Sequencing/TBISeq/addressing reviewer comments/mapped reference/azimuth_umap.Rds") #'azimuth_umap.Rds'
object <- refdata[, Cells(projected.umap)]
object[['umap.proj']] <- projected.umap
mappedref <- object

predictions <- read.delim("~/Documents/Sequencing/TBISeq/addressing reviewer comments/mapped reference/azimuth_pred.tsv", row.names = 1) #'azimuth_pred.tsv'
mappedref <- AddMetaData(
  object = mappedref,
  metadata = predictions)

ref_obj <- list(mappedref = mappedref)

ref_obj <- lapply(X = ref_obj, FUN = function(x) {
  x <- PercentageFeatureSet(object = x, pattern = "^mt-", col.name = "percent.mt")
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- NormalizeData(x)
  x <- SCTransform(object = x, vars.to.regress = c("percent.mt"), verbose = FALSE)
  x <- RunPCA(object = x, npcs = 50)
  x <- RunUMAP(object = x, dims = 1:30, verbose = FALSE, return.model = TRUE)
  x <- FindNeighbors(object = x, dims = 1:30, verbose = FALSE)
  x <- FindClusters(object = x, verbose = FALSE)
})

refdata_allen <- ref_obj[[1]]

DimPlot(refdata_allen, group.by = 'predicted.subclass', label = T)

```


```{r subset_refdata_neurons, echo=FALSE}

# Subset neurons from refdata
Idents(object = refdata_allen) <- refdata_allen@meta.data$'predicted.subclass'
refdata.neurons <- subset(refdata_allen, idents = c('Oligo', 'Astro', 'Micro-PVM', 'Endo'), invert = T)
refdata.neurons <- RunPCA(object = refdata.neurons,  ndims.print = 1:20, nfeatures.print = 12)
ElbowPlot(refdata.neurons, ndims=60)
refdata.neurons <- RunUMAP(object = refdata.neurons, reduction = "pca", dims = c(1:50))
refdata.neurons <- FindNeighbors(object = refdata.neurons, reduction = "pca", dims = c(1:50))
refdata.neurons <- FindClusters(refdata.neurons, resolution = 0.6)
refdata.neurons <- RunTSNE(object = refdata.neurons, reduction = "pca", dims = c(1:50))
DefaultAssay(object = refdata.neurons) <- "SCT"
DimPlot(refdata.neurons, group.by = 'predicted.subclass', label = TRUE, label.size = 4, repel = TRUE)

```


```{r refmap_atf3_lognorm, echo=FALSE}

DefaultAssay(atf3creer) <- 'integrated'
atf3creer <- RunUMAP(atf3creer, dims=1:30, return.model = T)

reference <- refdata_allen

DefaultAssay(atf3creer) <- 'RNA'
DefaultAssay(reference) <- 'RNA'

anchors <- FindTransferAnchors(
  reference = reference,
  query = atf3creer,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30
)

atf3creer_mapped <- MapQuery(
  anchorset = anchors,
  query = atf3creer,
  reference = reference,
  refdata = list(
    subclass_lognorm = "predicted.subclass",
    class_lognorm = "predicted.class"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

p2 = DimPlot(atf3creer_mapped, reduction = "ref.umap", label = TRUE, label.size = 4, repel = TRUE, group.by = 'predicted.subclass_lognorm') + NoLegend() + ggtitle('Atf3CreER') + ylim(-18,18) + NoAxes()
p3 = DimPlot(reference, group.by = 'predicted.subclass', label.size = 4, label = T, repel = T) + NoLegend() + ggtitle('Reference') + ylim(-18,18) + NoAxes()
plot_grid( p2, p3, ncol = 2)

atf3creer <- AddMetaData(object = atf3creer, metadata = atf3creer_mapped$predicted.subclass_lognorm, col.name = 'predicted.subclass_lognorm')

DimPlot(atf3creer, group.by = 'predicted.subclass_lognorm', label = TRUE, label.size = 4, repel = TRUE)

```


```{r refmap_atf3_sct, echo=FALSE}

reference <- refdata_allen

DefaultAssay(atf3creer) <- 'SCT'
DefaultAssay(reference) <- 'SCT'

anchors <- FindTransferAnchors(
  reference = reference,
  query = atf3creer,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30
)

atf3creer_mapped <- MapQuery(
  anchorset = anchors,
  query = atf3creer,
  reference = reference,
  refdata = list(
    subclass_sct = "predicted.subclass",
    class_sct = "predicted.class"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

p2 = DimPlot(atf3creer_mapped, reduction = "ref.umap", label = TRUE, label.size = 4, repel = TRUE, group.by = 'predicted.subclass_sct') + NoLegend() + ggtitle('Atf3CreER') + ylim(-18,18) + NoAxes()
p3 = DimPlot(reference, group.by = 'subclass', label.size = 4, label = T, repel = T) + NoLegend() + ggtitle('Reference') + ylim(-18,18) + NoAxes()
plot_grid( p2, p3, ncol = 2)

atf3creer <- AddMetaData(object = atf3creer, metadata = atf3creer_mapped$predicted.subclass_sct, col.name = 'predicted.subclass_sct')

DimPlot(atf3creer, group.by = 'predicted.subclass_sct', label = TRUE, label.size = 4, repel = TRUE)
DimPlot(atf3creer, group.by = 'predicted.subclass_lognorm', label = TRUE, label.size = 4, repel = TRUE)

## Assign weird group as "Unknown"
Idents(atf3creer) <- atf3creer@meta.data$"predicted.subclass_lognorm"
plot <- DimPlot(atf3creer)
atf3creer <- CellSelector(plot, object = atf3creer, ident = "Unknown")
atf3creer$predicted.subclass_lognorm_edit <- Idents(atf3creer)
DimPlot(atf3creer, group.by = 'predicted.subclass_lognorm_edit')

```


```{r subset_atf3_neurons, echo=FALSE}

neurons_Atf3CreER <- subset(atf3creer, idents = c('Oligo', 'Astro', 'Micro-PVM', 'VLMC', 'Endo', 'Peri', 'OPC', 'Unknown'), invert = T)

DefaultAssay(object = neurons_Atf3CreER) <- "integrated"
neurons_Atf3CreER <- RunPCA(object = neurons_Atf3CreER,  ndims.print = 1:20, nfeatures.print = 12)
ElbowPlot(neurons_Atf3CreER, ndims=60)
neurons_Atf3CreER <- RunUMAP(object = neurons_Atf3CreER, reduction = "pca", dims = c(1:40), n.neighbors = 35)
neurons_Atf3CreER <- FindNeighbors(object = neurons_Atf3CreER, dims = c(1:40))
neurons_Atf3CreER <- FindClusters(neurons_Atf3CreER, resolution = 0.6)
neurons_Atf3CreER <- RunTSNE(object = neurons_Atf3CreER, reduction = "pca", dims = c(1:40))
DefaultAssay(object = neurons_Atf3CreER) <- "SCT"


neurons_Atf3CreER$predicted.subclass_lognorm_edit <- factor(x = neurons_Atf3CreER$predicted.subclass_lognorm_edit,
                                                            levels = c('L2/3 IT', 'L5 IT', 'L5 ET', 'L5/6 NP', 'L6 IT', 'L6 CT', 'L6b',
                                                                       'Lamp5', 'Pvalb', 'Sncg', 'Sst', 'Vip'))
DimPlot(neurons_Atf3CreER, group.by = 'predicted.subclass_lognorm_edit')

```


```{r save_objects, echo=FALSE}

save(refdata, file = "Allen_mop_processed_21Jul24.Robj")
save(refdata.neurons, file = "Allen_mop_neurons_21Jul24.Robj")
save(atf3creer, file = "integrated_Atf3CreER_MF_mappedtoreference_23Jul24.Robj")
save(neurons_Atf3CreER, file = "integrated_neurons_Atf3CreER_MF_subsetatf3creer_umapn.neighbors_23Jul24.Robj")

```


```{r fig2plots, echo=FALSE}

DimPlot(neurons_Atf3CreER, group.by = 'predicted.subclass_lognorm_edit', pt.size = .1) + NoAxes()

Idents(neurons_Atf3CreER) <- neurons_Atf3CreER@meta.data$"predicted.subclass_lognorm_edit"
neurons_Atf3CreER_subset <- subset(neurons_Atf3CreER, idents = c('L5/6 NP', 'L6b', 'Sncg'), invert = T)

# save objects in a list
seurat.list <- list(
    reference=refdata.neurons,
    atf3=neurons_Atf3CreER_subset)

# Assign column names for annotation info in each metadata
annocol.list <- list(
    reference='predicted.subclass',
    atf3='predicted.subclass_lognorm_edit'
)

# Add more info to each metadata
seurat.list <- map(names(seurat.list), function(name) {
    obj <- seurat.list[[name]]
    obj$group <- name  # adds dataset info to metadata 
    obj$type_group <- obj@meta.data[[annocol.list[[name]]]]  # adds celltype annotation to metadata
    return(obj)
}) %>%
set_names(names(seurat.list))

# Prep a vector for common gene symbols
common.genes <- intersect(
    rownames(seurat.list[[1]]),
    rownames(seurat.list[[2]])
)

# Prep a vector for common column names in the metadata
common.columns <- intersect(
    colnames(seurat.list[[1]]@meta.data),
    colnames(seurat.list[[2]]@meta.data)
)

# Create a list for objects having subsetted metadata and genes
subset.list <- lapply(seurat.list, function(obj) {
    # Subset metadata by column
    obj@meta.data <- obj@meta.data[, common.columns]
    # Subset SCT counts by gene
    mat <- subset(obj[['SCT']], features=common.genes)
    # Build a new Seurat obj using the subsetted counts and metadata
    obj_new <- CreateSeuratObject(
        counts=mat,
        meta.data=obj@meta.data
    )
    return(obj_new)
})

subset.obj <- merge(subset.list[[1]], subset.list[[2]])

###### Issues with subscript out of bounds error - Mira ran same code above and it worked and sent me the object
merged_obj <- readRDS("../addressing reviewer comments/seurat_merged.rds")

Idents(merged_obj) <- merged_obj@meta.data$"type_group"
merged_obj_subset <- subset(merged_obj, idents = c('OPC', 'Peri', 'VLMC', 'L5/6 NP', 'L6 IT Car3',
                                                   'L6b', 'Sst Chodl', 'Sncg'), invert = T)
merged_obj_subset$group_type <- paste0(merged_obj_subset$group, '_', merged_obj_subset$type_group)

# refdata.neurons$group <- 'reference'
# refdata.neurons$type_group <- paste0(refdata.neurons$predicted.subclass, '_', refdata.neurons$group)
# neurons_Atf3CreER_subset$type_group <- paste0(neurons_Atf3CreER_subset$predicted.subclass_lognorm_edit, '_', neurons_Atf3CreER$group)
# 
# genes <- rownames(neurons_Atf3CreER_subset@assays$SCT@data)
# subset_refdata.neurons <- subset(refdata.neurons, features = genes)
# genes <- rownames(subset_refdata.neurons@assays$SCT@data)
# subset_atf3creer_neurons <- subset(neurons_Atf3CreER_subset, features = genes)
# 
# merge <- merge(subset_refdata.neurons, neurons_Atf3CreER_subset)
Idents(merged_obj_subset) <- 'group_type'
# DefaultAssay(merged_obj) <- 'SCT'
# DotPlot(merge, features = corticalneurons, group.by = 'type_group') + RotatedAxis() + ggtitle("SCT assay")

corticalneurons <- c('Slc17a7',
                     'Cux2', 'Rorb',
                     'Fezf2', 'Sulf1', 'Bcl11b',
                     'Foxp2',
                     'Gad1', 'Gad2',
                     'Lamp5', 'Pvalb', 'Sst', 'Vip')
DotPlot(merged_obj_subset, features = corticalneurons, group.by = 'group_type') + RotatedAxis()

stressgenes <- c('Ecel1', 'Sprr1a', 'Sox11', 'Ret', 'Nmnat2', 'Atf3',
                 'Rest', 'Creb5', 'Jun', 'Oxsr1', 'Stat3',
                 'Bbc3', 'Casp3', 'Htt', 'Ddit3',
                 'Eif2ak3', 'Ern1')
DotPlot(neurons_Atf3CreER_subset, features = stressgenes, group.by = 'predicted.subclass_lognorm_edit') + RotatedAxis()


```

```{r fig4plots, echo=FALSE}

ionchannels <- c('Scn1a', 'Scn2a', 'Scn3a',
                 'Kcnq2', 'Kcnq3', 'Kcnq5',
                 'Kcnb1', 'Kcnab1', 'Cacna1a', 'Hcn1',
                 'Atp1a1', 'Atp1b1', 'Kcnma1', 'Kcnk2')
Idents(merged_obj_subset) <- 'type_group'
merged_obj_L2and5 <- subset(merged_obj_subset, idents = c('L2/3 IT', 'L5 ET', 'L5 IT'))
DotPlot(merged_obj_L2and5, features = ionchannels, group.by = 'group_type') + RotatedAxis() + coord_flip()


```

```{r suppfig2plots, echo=FALSE}

atf3creer$predicted.subclass_lognorm_edit <- factor(x = atf3creer$predicted.subclass_lognorm_edit,
                                                            levels = c('Astro', 'Endo', 'Micro-PVM', 'Oligo', 'VLMC', 'Peri', 'OPC',
                                                                       'L2/3 IT', 'L5 IT', 'L5 ET', 'L5/6 NP', 'L6 IT',
                                                                       'L6 CT', 'L6b','Lamp5', 'Pvalb', 'Sncg', 'Sst', 'Vip',
                                                                       'Unknown'))
DimPlot(atf3creer, group.by = 'predicted.subclass_lognorm_edit')

atf3creer_subset <- subset(atf3creer, idents = c('Unknown'), invert = T)
Idents(object = atf3creer_subset) <- atf3creer_subset@meta.data$'predicted.subclass_lognorm_edit'
atf3creer_subset$subtypes <- Idents(atf3creer_subset)
current_ids <- unique(atf3creer_subset@meta.data$predicted.subclass_lognorm_edit)
atf3creer_subset@meta.data$subtypes <- plyr::mapvalues(x = atf3creer_subset@meta.data$subtypes,
                                                            from = current_ids,
                                                            to = c('Neurons', 'Micro-PVM', 'Neurons', 'Neurons', 'Neurons',
                                                                   'VLMC', 'Neurons', 'Neurons', 'Astro', 'Oligo', 'Neurons',
                                                                   'Endo', 'Peri', 'OPC', 'Neurons', 'Neurons', 'Neurons',
                                                                   'Neurons', 'Neurons'))                                                                    

DimPlot(atf3creer_subset, group.by = 'subtypes')

atf3creer_subset$subtypes <- factor(x = atf3creer_subset$subtypes, levels = c('Neurons', 'Micro-PVM', 'Astro', 'VLMC',
                                                                              'Peri', 'Endo', 'Oligo', 'OPC'))

corticalcells <- c('Snap25', 'Syt7',
                   'Ctss', 'Cx3cr1', 'C1qb',
                   'Gfap', 'Aqp4',
                   'Pdgfra',
                   'Bgn',
                   'Sulf1',
                   'Mbp', 'Mag', 'Mog', 'Cd9')
DotPlot(atf3creer_subset, features = corticalcells, group.by = 'subtypes') + RotatedAxis()

```


```{r suppfig3plots, echo=FALSE}

FeaturePlot(neurons_Atf3CreER, c('Atf3', 'Ecel1', 'Ddit3', 'Satb2', 'Gad2'), ncol = 5)

```


