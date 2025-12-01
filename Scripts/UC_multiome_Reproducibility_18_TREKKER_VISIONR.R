#####################################################
# TREKKER VISIONR
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

Get_signature_score_literature <- function(lineage){
  if(lineage == 'TAM'){
    DOGMA = readRDS(file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_Myeloid.rds")))
    Idents(DOGMA) = 'coarse_celltype'
    DOGMA = subset(DOGMA, idents = c('Mono_Mac'))
  }else if(lineage == 'BCG'){
    DOGMA = readRDS(file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_BCG.rds")))
    Idents(DOGMA) = 'coarse_celltype'
    DOGMA = subset(DOGMA, idents = c('Mono_Mac'))
  }else{
    DOGMA = readRDS(file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")))
  }
  DefaultAssay(DOGMA) = 'RNA'
  gene_df = fread(paste0("Reproducibility/Data/Signature_gene_list/TableS3_module_genes_literature_",lineage,".txt")) %>% as.data.frame()

  # Convert data.frame / data.table to list of columns
  sig_list <- lapply(gene_df, function(x) {
    x <- as.character(x)
    x <- x[x != ""]        # remove empty strings
    x <- na.omit(x)        # remove NA
    unique(x)
  })

  signatures=c()
  for(tmp_num in 1:length(sig_list)){
      tmp_genes = sig_list[[tmp_num]]
      tmp_genes2 = intersect(tmp_genes, rownames(DOGMA))
      tmp_module = c(rep(1, length(tmp_genes2)))
      names(tmp_module) = tmp_genes2
      sig <- createGeneSignature(name = names(sig_list)[[tmp_num]], sigData = tmp_module)
      signatures = c(signatures,sig)
  }
  
  vision.obj <- Vision(DOGMA,
                       signatures = signatures,
                       dimRed = "scvi_latent",
                       projection_methods = c("UMAP"),
                       meta = DOGMA@meta.data)
  
  options(mc.cores = 4)
  vis <- analyze(vision.obj)
  
  sig_df = getSignatureScores(vis) %>% as.data.frame()
  write.table_n_2(sig_df, 'CB', paste0('Reproducibility/Results/VISIONR/UC_DOGMA_',lineage,'_signature_score_literature.txt'))
  }

Get_signature_score_hotspot <- function(lineage){
  DOGMA = readRDS(file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")))

  DefaultAssay(DOGMA) = 'RNA'
  gene_df = fread(paste0("Reproducibility/Data/Signature_gene_list/TableS3_module_genes_",lineage,"_hotspot.txt")) %>% as.data.frame()
  gene_df <- gene_df[gene_df$Annotation != "" & !is.na(gene_df$Annotation), ]

  sig_list <- gene_df %>%
  group_by(Annotation) %>%
  summarise(genes = list(Gene)) %>%
  deframe()

  names(sig_list) <- gsub(" ", "_", names(sig_list))

  signatures=c()
  for(tmp_num in 1:length(sig_list)){
      tmp_genes = sig_list[[tmp_num]]
      tmp_genes2 = intersect(tmp_genes, rownames(DOGMA))
      tmp_module = c(rep(1, length(tmp_genes2)))
      names(tmp_module) = tmp_genes2
      sig <- createGeneSignature(name = names(sig_list)[[tmp_num]], sigData = tmp_module)
      signatures = c(signatures,sig)
  }
  
  vision.obj <- Vision(DOGMA,
                       signatures = signatures,
                       dimRed = "scvi_latent",
                       projection_methods = c("UMAP"),
                       meta = DOGMA@meta.data)
  
  options(mc.cores = 4)
  vis <- analyze(vision.obj)
  
  sig_df = getSignatureScores(vis) %>% as.data.frame()
  write.table_n_2(sig_df, 'CB', paste0('Reproducibility/Results/VISIONR/UC_DOGMA_',lineage,'_signature_score_hotspot.txt'))
  }

Get_signature_score_scenic_module <- function(lineage){
  DOGMA = readRDS(file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_CD8_T_NK_ILC.rds")))
  DefaultAssay(DOGMA) = 'RNA'
  gene_df = fread(paste0("Reproducibility/Data/Signature_gene_list/TableS3_module_genes_CD8_T_scenicplus.txt")) %>% as.data.frame()

  # Convert data.frame / data.table to list of columns
  sig_list <- lapply(gene_df, function(x) {
    x <- as.character(x)
    x <- x[x != ""]        # remove empty strings
    x <- na.omit(x)        # remove NA
    unique(x)
  })

  signatures=c()
  for(tmp_num in 1:length(sig_list)){
      tmp_genes = sig_list[[tmp_num]]
      tmp_genes2 = intersect(tmp_genes, rownames(DOGMA))
      tmp_module = c(rep(1, length(tmp_genes2)))
      names(tmp_module) = tmp_genes2
      sig <- createGeneSignature(name = names(sig_list)[[tmp_num]], sigData = tmp_module)
      signatures = c(signatures,sig)
  }
  
  vision.obj <- Vision(DOGMA,
                       signatures = signatures,
                       dimRed = "scvi_latent",
                       projection_methods = c("UMAP"),
                       meta = DOGMA@meta.data)
  
  options(mc.cores = 4)
  vis <- analyze(vision.obj)
  
  sig_df = getSignatureScores(vis) %>% as.data.frame()
  write.table_n_2(sig_df, 'CB', paste0('Reproducibility/Results/VISIONR/UC_DOGMA_CD8_T_NK_ILC_signature_score_scenicplue_module.txt'))
  }

Get_signature_score_literature(lineage = 'TAM')
Get_signature_score_literature(lineage = 'CD4_T')
Get_signature_score_literature(lineage = 'MSC')
Get_signature_score_literature(lineage = 'Malignant')
Get_signature_score_literature(lineage = 'BCG')
Get_signature_score_hotspot(lineage = 'MSC')
Get_signature_score_hotspot(lineage = 'Malignant')
Get_signature_score_scenic_module(lineage = 'CD8_T_NK_ILC')

