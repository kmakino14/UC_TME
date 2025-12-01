#####################################################
# 03. DOGMA Seurat object construction
#####################################################

setwd(path_to_wd)
source("Reproducibility/Scripts/Source/my.source.R")
source("Reproducibility/Scripts/Source/Seurat_source.R")
options(stringsAsFactors=F)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(cicero))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(VISION))

suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(ggrepel))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsci))
suppressMessages(library(colorspace))
suppressMessages(library(BuenColors))

suppressMessages(library(survival))
suppressMessages(library(survminer))

suppressMessages(library(dplyr))
set.seed(1234)

#******************************************
data_dir="Reproducibility/Data"

# Function to build Seurat object for a given lineage
build_seurat_object <- function(lineage, data_dir = "Reproducibility/Data") {
  
  message("Processing lineage: ", lineage)
  
  # Input files
  frag.file   <- file.path(data_dir, "ATAC_fragments", "atac_fragments_Global.tsv.gz")
  exp_matrix  <- Read10X(file.path(data_dir, "matrix_files", lineage, "GEX"))
  ATAC_counts <- Read10X(file.path(data_dir, "matrix_files", lineage, "SnapATAC2"))
  
  # Filter ATAC counts to standard chromosomes
  grange.counts <- StringToGRanges(rownames(ATAC_counts), sep = c(":", "-"))
  grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  ATAC_counts   <- ATAC_counts[as.vector(grange.use), ]
  
  # Gene annotations
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  
  # Chromatin assay
  chrom_assay <- CreateChromatinAssay(
    counts     = ATAC_counts,
    sep        = c(":", "-"),
    genome     = "hg38",
    fragments  = frag.file,
    min.cells  = 10,
    annotation = annotations
  )
  
  # Metadata
  meta_data_all <- fread_n(file.path(data_dir, "UC_DOGMA_metadata.txt"))
  common_cells  <- intersect(rownames(meta_data_all), colnames(chrom_assay))
  meta_data     <- meta_data_all[common_cells, ]
  
  # Seurat object
  DOGMA <- CreateSeuratObject(
    counts    = subset(chrom_assay, cells = common_cells),
    assay     = "ATAC",
    meta.data = meta_data
  )
  
  # Add RNA assay
  DOGMA[["RNA"]] <- CreateAssayObject(counts = exp_matrix)
  
  # RNA preprocessing
  DefaultAssay(DOGMA) <- "RNA"
  DOGMA <- NormalizeData(DOGMA)
  DOGMA <- FindVariableFeatures(DOGMA, selection.method = "vst", nfeatures = 3000)
  DOGMA <- ScaleData(DOGMA, features = rownames(DOGMA))
  DOGMA <- RunPCA(DOGMA, features = VariableFeatures(DOGMA))
  
  # ATAC preprocessing
  DefaultAssay(DOGMA) <- "ATAC"
  DOGMA <- FindTopFeatures(DOGMA, min.cutoff = 5)
  DOGMA <- RunTFIDF(DOGMA)
  DOGMA <- RunSVD(DOGMA)
  
  # Add ADT assay
  adt_counts_all <- fread_n(file.path(data_dir, "UC_DOGMA_ADT_counts.txt")) %>% t() %>% as.matrix()
  adt_counts <- adt_counts_all[, colnames(DOGMA), drop = FALSE]
  DOGMA[["ADT"]] <- CreateAssayObject(counts = adt_counts)
  DefaultAssay(DOGMA) <- "ADT"
  DOGMA <- NormalizeData(DOGMA, normalization.method = "CLR", margin = 2)
  DOGMA <- ScaleData(DOGMA, assay = "ADT")
  
  # Add embeddings
  latent_df <- fread_n(file.path(data_dir, "embeddings", paste0("UC_DOGMA_", lineage, "_latent_df.txt")))
  DOGMA[["scvi_latent"]] <- CreateDimReducObject(embeddings = as.matrix(latent_df), key = "scvi_", assay = "ATAC")
  
  umap_df <- fread_n(file.path(data_dir, "embeddings", paste0("UC_DOGMA_", lineage, "_UMAP_df.txt")))
  DOGMA[["UMAP"]] <- CreateDimReducObject(embeddings = as.matrix(umap_df), key = "umap_", assay = "ATAC")
  
  # Save RDS
  saveRDS(DOGMA, file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")))
}

build_seurat_object(lineage = 'CD4_T',        data_dir = data_dir)
build_seurat_object(lineage = 'CD8_T_NK_ILC', data_dir = data_dir)
build_seurat_object(lineage = 'B',            data_dir = data_dir)
build_seurat_object(lineage = 'Myeloid',      data_dir = data_dir)
build_seurat_object(lineage = 'MSC',          data_dir = data_dir)
build_seurat_object(lineage = 'Malignant',    data_dir = data_dir)
build_seurat_object(lineage = 'BCG',          data_dir = data_dir)