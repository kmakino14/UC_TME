#####################################################
# DOGMA Cicero coaccesibility
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
suppressMessages(library(dplyr))
set.seed(1234)

##################################
data_dir="Reproducibility/Data"

lineage_list = c('Malignant', 'CD4_T', 'CD8_T_NK_ILC', 'B', 'Myeloid', 'MSC' ,'BCG')
lineage = 'Malignant' ## Choose one from lineage_list

DOGMA = readRDS(file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")))
	
DefaultAssay(DOGMA) = 'ATAC'
DOGMA.cds <- as.CellDataSet(DOGMA)
umap_df <- Embeddings(DOGMA, "UMAP") %>% as.data.frame()
DOGMA.cicero <- make_cicero_cds(DOGMA.cds, reduced_coordinates = umap_df)
	
# get the chromosome sizes from the Seurat object
genome_all <- seqlengths(DOGMA)

for(k in 1:23){
# convert chromosome sizes to a dataframe
genome <- genome_all[k]
tmp_chr <- names(genome)
genome.df <- data.frame("chr" = tmp_chr, "length" = genome)

# run cicero
conns <- run_cicero(DOGMA.cicero, genomic_coords = genome.df, sample_num = 100)
ccans <- generate_ccans(conns, coaccess_cutoff_override = .1)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
saveRDS(links, file = paste0("Reproducibility/Results/Cicero/",lineage,"/UC_DOGMA_links_",lineage,"_",tmp_chr,".rds"))
rm(conns, ccans, links)
gc(verbose = FALSE)
}