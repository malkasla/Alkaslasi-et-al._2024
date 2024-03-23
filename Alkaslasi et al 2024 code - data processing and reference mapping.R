###### Alkaslasi et al. 2024 ######

#### Female data pre-processing ####

# Import data
Atf3CreER_7dpi_F.data <- Read10X(data.dir = "Atf3CreER")

# Save un-processed data
save(Atf3CreER_7dpi_F.data, file = "Atf3CreER_7dpi_F.data.robj")

colnames(Atf3CreER_7dpi_F.data) <- paste0(colnames(Atf3CreER_7dpi_F.data),"-Atf3CreER_TBI")
Atf3CreER_7dpi_F <-CreateSeuratObject(Atf3CreER_7dpi_F.data)
Atf3CreER_7dpi_F <- PercentageFeatureSet(object = Atf3CreER_7dpi_F, pattern = "^mt-", col.name = "percent.mt")
Atf3CreER_7dpi_F <- SCTransform(object = Atf3CreER_7dpi_F, vars.to.regress = "percent.mt", verbose = FALSE)
Atf3CreER_7dpi_F <- RunPCA(object = Atf3CreER_7dpi_F, npcs = 50)
ElbowPlot(Atf3CreER_7dpi_F, ndims = 50)
Atf3CreER_7dpi_F <- RunUMAP(object = Atf3CreER_7dpi_F, dims = 1:35, verbose = FALSE)
Atf3CreER_7dpi_F <- FindNeighbors(object = Atf3CreER_7dpi_F, dims = 1:35, verbose = FALSE)
Atf3CreER_7dpi_F <- FindClusters(object = Atf3CreER_7dpi_F, verbose = FALSE, resolution = 0.6)
DimPlot(object = Atf3CreER_7dpi_F, label = TRUE) + NoLegend()

Idents(object = Atf3CreER_7dpi_F) <- Atf3CreER_7dpi_F@meta.data$orig.ident
VlnPlot(Atf3CreER_7dpi_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Atf3CreER_7dpi_F, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Atf3CreER_7dpi_F, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
Atf3CreER_7dpi_F <- subset(Atf3CreER_7dpi_F, subset = percent.mt < 10)

Atf3CreER_7dpi_F$sex <- 'F'
save(Atf3CreER_7dpi_F, file = "Atf3CreER_7dpi_F.robj")

#### Male data pre-processing ####

# Import data - filtered files
Atf3CreER_7dpi_M.data <- Read10X(data.dir = "Atf3CreER_M")

colnames(Atf3CreER_7dpi_M.data) <- paste0(colnames(Atf3CreER_7dpi_M.data),"-Atf3CreER_TBI")
Atf3CreER_7dpi_M <-CreateSeuratObject(Atf3CreER_7dpi_M.data)
Atf3CreER_7dpi_M <- PercentageFeatureSet(object = Atf3CreER_7dpi_M, pattern = "^mt-", col.name = "percent.mt")

Idents(object = Atf3CreER_7dpi_M) <- Atf3CreER_7dpi_M@meta.data$orig.ident
VlnPlot(Atf3CreER_7dpi_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Atf3CreER_7dpi_M, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Atf3CreER_7dpi_M, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
Atf3CreER_7dpi_M <- subset(Atf3CreER_7dpi_M, subset = percent.mt < 10 & nCount_RNA < 20000)

Atf3CreER_7dpi_M <- SCTransform(object = Atf3CreER_7dpi_M, vars.to.regress = "percent.mt", verbose = FALSE)
Atf3CreER_7dpi_M <- RunPCA(object = Atf3CreER_7dpi_M, npcs = 50)
ElbowPlot(Atf3CreER_7dpi_M, ndims = 50)
Atf3CreER_7dpi_M <- RunUMAP(object = Atf3CreER_7dpi_M, dims = 1:35, verbose = FALSE)
Atf3CreER_7dpi_M <- FindNeighbors(object = Atf3CreER_7dpi_M, dims = 1:35, verbose = FALSE)
Atf3CreER_7dpi_M <- FindClusters(object = Atf3CreER_7dpi_M, verbose = FALSE, resolution = 0.6)
DimPlot(object = Atf3CreER_7dpi_M, label = TRUE) + NoLegend()

Atf3CreER_7dpi_M$sex <- 'M'

#### Integrate datasets ####

control.anchors <-FindIntegrationAnchors(object.list = list(Atf3CreER_7dpi_F, Atf3CreER_7dpi_M), dims = 1:30)
Atf3CreER_int_MF_7dpi <- IntegrateData(anchorset = control.anchors, dims = 1:30)
DefaultAssay(object = Atf3CreER_int_MF_7dpi) <- "integrated"
Atf3CreER_int_MF_7dpi <- ScaleData(object = Atf3CreER_int_MF_7dpi, verbose = FALSE)
Atf3CreER_int_MF_7dpi <- RunPCA(object = Atf3CreER_int_MF_7dpi, npcs = 30, verbose = FALSE)
ElbowPlot(Atf3CreER_int_MF_7dpi, ndims = 40)
Atf3CreER_int_MF_7dpi <- RunUMAP(object = Atf3CreER_int_MF_7dpi, reduction = "pca", dims = 1:30)
Atf3CreER_int_MF_7dpi <- FindNeighbors(object = Atf3CreER_int_MF_7dpi, reduction = "pca", dims = 1:30)
Atf3CreER_int_MF_7dpi <- FindClusters(Atf3CreER_int_MF_7dpi, resolution = 0.8)
Atf3CreER_int_MF_7dpi <- RunTSNE(object = Atf3CreER_int_MF_7dpi, reduction = "pca", dims = c(1:30))
DefaultAssay(object = Atf3CreER_int_MF_7dpi) <- "SCT"

#### Reference mapping ####

# Load reference data - AllenSampleSubset_Mop - Azimuth
reference <- refdata

DefaultAssay(Atf3CreER_int_MF_7dpi) <- 'SCT'
DefaultAssay(refdata) <- 'SCT'

anchors <- FindTransferAnchors(
  reference = reference,
  query = Atf3CreER_int_MF_7dpi,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30)

Atf3CreER_int_MF_7dpi_mapped <- MapQuery(
  anchorset = anchors,
  query = Atf3CreER_int_MF_7dpi,
  reference = reference,
  refdata = list(
    subclass = "subclass",
    class = "class"),
  reference.reduction = "pca", 
  reduction.model = "umap")

Atf3CreER_int_MF_7dpi <- AddMetaData(object = Atf3CreER_int_MF_7dpi, metadata = Atf3CreER_int_MF_7dpi_mapped$predicted.subclass, col.name = 'predicted.subclass')

# Define unknown cluster - expresses hippocampal markers
Idents(Atf3CreER_int_MF_7dpi) <- Atf3CreER_int_MF_7dpi@meta.data$"predicted.subclass"
plot <- DimPlot(Atf3CreER_int_MF_7dpi)
Atf3CreER_int_MF_7dpi <- CellSelector(plot, object = Atf3CreER_int_MF_7dpi, ident = "Unknown")
Atf3CreER_int_MF_7dpi$predicted.subclass_edit <- Idents(Atf3CreER_int_MF_7dpi)
DimPlot(Atf3CreER_int_MF_7dpi, group.by = 'predicted.subclass_edit')

#### Subset neurons ####

# From reference data
refdata.neurons <- subset(refdata, idents = c('21', '22', '25', '19', '5', '20', '26', '17'), invert = T)
refdata.neurons <- RunPCA(object = refdata.neurons,  ndims.print = 1:20, nfeatures.print = 12)
ElbowPlot(refdata.neurons, ndims=60)
refdata.neurons <- RunUMAP(object = refdata.neurons, reduction = "pca", dims = c(1:40), return.model = T)
refdata.neurons <- FindNeighbors(object = refdata.neurons, reduction = "pca", dims = c(1:40))
refdata.neurons <- FindClusters(refdata.neurons, resolution = 0.6)
refdata.neurons <- RunTSNE(object = refdata.neurons, reduction = "pca", dims = c(1:40))

# From integrated data
Atf3CreER_neurons_MF <- subset(Atf3CreER_int_MF_7dpi, idents = c('6', '8', '19', '12', '3', '17', '24', '20', '14', '4'))
DefaultAssay(object = Atf3CreER_neurons_MF) <- "integrated"
Atf3CreER_neurons_MF <- RunPCA(object = Atf3CreER_neurons_MF, npcs = 30)
ElbowPlot(Atf3CreER_neurons_MF,ndims=40)
Atf3CreER_neurons_MF <- RunUMAP(object = Atf3CreER_neurons_MF, reduction = "pca", dims = c(1:20))
Atf3CreER_neurons_MF <- FindNeighbors(object = Atf3CreER_neurons_MF, reduction = "pca", dims = c(1:20))
Atf3CreER_neurons_MF <- FindClusters(Atf3CreER_neurons_MF, resolution = 0.6) ## Mess with resolution
Atf3CreER_neurons_MF <- RunTSNE(object = Atf3CreER_neurons_MF, reduction = "pca", dims = c(1:20))
DefaultAssay(object = Atf3CreER_neurons_MF) <- "SCT"

#### Reference map neurons ####

reference <- refdata.neurons

DefaultAssay(Atf3CreER_subset_neurons_MF) <- 'SCT'
DefaultAssay(reference) <- 'SCT'
anchors <- FindTransferAnchors(
  reference = reference,
  query = Atf3CreER_subset_neurons_MF,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30)

Atf3CreER_subset_neurons_MF_mapped <- MapQuery(
  anchorset = anchors,
  query = Atf3CreER_subset_neurons_MF,
  reference = reference,
  refdata = list(
    subclass = "subclass",
    class = "class"),
  reference.reduction = "pca", 
  reduction.model = "umap")

Atf3CreER_subset_neurons_MF <- AddMetaData(object = Atf3CreER_subset_neurons_MF, metadata = Atf3CreER_subset_neurons_MF_mapped$predicted.subclass, col.name = 'predicted.subclass')



