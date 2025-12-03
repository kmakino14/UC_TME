#####################################################
# DOGMA DEG analysis by MAST
#####################################################

setwd(path_to_wd)
source("Reproducibility/Scripts/Source/my.source.R")
source("Reproducibility/Scripts/Source/Seurat_source.R")
options(stringsAsFactors=F)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(MAST))
set.seed(1234)

#******************************************
data_dir="Reproducibility/Data"

###################################
## DEG each celltype
###################################
DEG_MAST <- function(lineage) {
  DOGMA = paste0("Reproducibility/Data/Seurat/UC_DOGMA_seurat_obj_",lineage,".rds") %>% readRDS()
  
  DefaultAssay(DOGMA) = "RNA"
  scRNA = DietSeurat(DOGMA, assays = 'RNA')
  
  Idents(scRNA) = 'RNA_HQ'
  scRNA = subset(scRNA, idents = c('HQ'))
  scRNA <- NormalizeData(scRNA)
  
  Idents(scRNA) = "celltype"
  deg <- Seurat::FindAllMarkers(object = scRNA,
                                only.pos = T, 
                                logfc.threshold = 0.5,
                                min.pct = 0,
                                return.thresh = 0.05,
                                slot = "data",
                                test.use = "MAST"
                                )
  rownames(deg) <- NULL
  deg = dplyr::select(deg, c("gene","avg_log2FC","p_val_adj","cluster"))
  write.table_n_2(deg, "No", paste0("Reproducibility/Results/Differential/UC_DOGMA_DEG_MAST_",lineage,"_v5.txt")) ### v5!
}

DEG_MAST(lineage = 'Malignant')
DEG_MAST(lineage = 'MSC')
DEG_MAST(lineage = 'CD4_T')
DEG_MAST(lineage = 'CD8_T_NK_ILC')
DEG_MAST(lineage = 'B')
DEG_MAST(lineage = 'Myeloid')


###################################
## DEG between LUM and NRP
###################################
DOGMA = paste0("Reproducibility/Data/Seurat/UC_DOGMA_seurat_obj_Malignant.rds") %>% readRDS()

DefaultAssay(DOGMA) = "RNA"
scRNA = DietSeurat(DOGMA, assays = 'RNA')
  
Idents(scRNA) = 'RNA_HQ'
scRNA = subset(scRNA, idents = c('HQ'))
scRNA <- NormalizeData(scRNA)

Idents(scRNA) = 'celltype'
Malignant = subset(scRNA, idents = c('NRP','LUM'))
Malignant <- NormalizeData(Malignant)

Idents(Malignant) = "celltype"
deg <- Seurat::FindMarkers(object = Malignant,
                           only.pos = F, 
                           logfc.threshold = 0,
                           min.pct = 0,
                           ident.1 = 'NRP',
                           ident.2 = 'LUM',
                           slot = "data",
                           test.use = "MAST"
                           )

deg = dplyr::select(deg, c("avg_log2FC","p_val_adj"))
write.table_n_2(deg, "No", paste0("Reproducibility/Results/Differential/UC_DOGMA_DEG_MAST_Malignant_NRP_vs_LUM.txt"))


###################################
## DEG between eTregs by STAGE
###################################
DOGMA = paste0("Reproducibility/Data/Seurat/UC_DOGMA_seurat_obj_CD4_T.rds") %>% readRDS()

DefaultAssay(DOGMA) = "RNA"
scRNA = DietSeurat(DOGMA, assays = 'RNA')
  
Idents(scRNA) = 'RNA_HQ'
scRNA = subset(scRNA, idents = ('HQ'))
scRNA <- NormalizeData(scRNA)

###################################
# RNA
metadata = fread_n('Reproducibility/Results/Milo/output/Milo_Treg_effector_metadata.txt')

Idents(scRNA) = "celltype"
CD4 = subset(scRNA, cells = rownames(metadata))
CD4 <- NormalizeData(CD4)
CD4$nhood_groups = metadata[colnames(CD4),]$nhood_groups

Idents(CD4) = "nhood_groups"
deg <- Seurat::FindMarkers(object = CD4,
                           only.pos = F, 
                           logfc.threshold = 0,
                           min.pct = 0,
                           ident.1 = "in_nhoods_MI_enr",
                           ident.2 = "in_nhoods_NMI_enr",
                           slot = "data",
                           test.use = "MAST"
                           )

write.table_n_2(deg, "Gene", paste0("Reproducibility/Results/Differential/UC_DOGMA_DEG_MAST_Treg_effector.txt"))

###################################
# ADT
counts = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_ADT_counts.txt") %>% t() %>% as.matrix()
CD4[['ADT']] = CreateAssayObject(counts = counts[,colnames(CD4)])
CD4 <- NormalizeData(CD4, normalization.method = "CLR", margin = 2)

DefaultAssay(CD4) = "ADT"
Idents(CD4) = "nhood_groups"
deg <- Seurat::FindMarkers(object = CD4,
                           only.pos = F, 
                           logfc.threshold = 0,
                           min.pct = 0,
                           ident.1 = "in_nhoods_MI_enr",
                           ident.2 = "in_nhoods_NMI_enr",
                           slot = "data",
                           test.use = "MAST"
                           )

write.table_n_2(deg, "Gene", paste0("Reproducibility/Results/Differential/UC_DOGMA_DEG_MAST_Treg_effector_ADT.txt"))

###################################
# TF
TF_zscore = fread_n('Reproducibility/Results/LINGER/Primary/cell_population_TF_activity_zscore_CD4_T.txt') %>% t()
CD4_TF = CreateSeuratObject(count = TF_zscore, data = TF_zscore, assay = "TF")

metadata = fread_n('Reproducibility/Results/Milo/output/Milo_Treg_effector_metadata.txt')
Treg = subset(CD4_TF, cells = rownames(metadata))
Treg$nhood_groups = metadata[colnames(Treg),]$nhood_groups

Idents(Treg) = "nhood_groups"
deg <- Seurat::FindMarkers(object = Treg,
                           only.pos = F, 
                           logfc.threshold = 0,
                           min.pct = 0,
                           ident.1 = "in_nhoods_MI_enr",
                           ident.2 = "in_nhoods_NMI_enr",
                           slot = "data",
                           test.use = "wilcox"
                           )
write.table_n_2(deg, "Gene", paste0("Reproducibility/Results/Differential/UC_DOGMA_DEG_wilcox_Treg_effector_TF.txt"))

###################################
## DEG between eTregs by State
###################################
df <- read.csv("Reproducibility/Results/CellRank2/CD4_T/fate_probabilities_dpt_Treg.csv", row.names = 1)
prob_score = df[, c("Treg_effector_2", "Treg_effector_3")] 

# Get barcodes
top_barcodes <- dplyr::arrange(prob_score, desc(Treg_effector_2)) %>% head(2000) %>% rownames()
bottom_barcodes <- dplyr::arrange(prob_score, desc(Treg_effector_3)) %>% head(2000) %>% rownames()

###################################
DOGMA = paste0("Reproducibility/Data/Seurat/UC_DOGMA_seurat_obj_CD4_T.rds") %>% readRDS()

DefaultAssay(DOGMA) = "RNA"
DOGMA = DietSeurat(DOGMA, assays = c('RNA','ADT'))

Idents(DOGMA) = "celltype"
Treg = subset(DOGMA, ident= c('Treg_naive', 'Treg_effector'))
Treg$state <- ifelse(
  colnames(Treg) %in% top_barcodes, "State1",
  ifelse(colnames(Treg) %in% bottom_barcodes, "State2", 'Other')
)

###################################
Idents(Treg) = 'RNA_HQ'
scRNA = subset(Treg, idents = ('HQ'))
scRNA <- NormalizeData(scRNA)

Idents(scRNA) = "state"
deg <- Seurat::FindMarkers(object = scRNA,
                           only.pos = F, 
                           logfc.threshold = 0,
                           min.pct = 0,
                           ident.1 = "State1",
                           ident.2 = "State2",
                           slot = "data",
                           test.use = "MAST"
                           )

write.table_n_2(deg, "Gene", paste0("Reproducibility/Results/Differential/UC_DOGMA_DEG_MAST_Treg_by_pseudotime_state.txt"))


###################################
## DEG between pre- and post-BCG
###################################
scRNA = readRDS("Reproducibility/Data/Seurat/UC_DOGMA_seurat_obj_Global_RNA_HQ.rds")

Idents(scRNA) = "sample"
BCG = subset(scRNA, idents = c("BC_011","BC_016","BC_023","BC_033","BC_037",
                               "BC_027","BC_032","BC_039","BC_040","BC_043","BC_044","BC_047","BC_048"))

BCG$pre_pos = fct_collapse(BCG$sample,
                           pre = c("BC_011","BC_016","BC_023","BC_033","BC_037"),
                           post = c("BC_027","BC_032","BC_039","BC_040","BC_043","BC_044","BC_047","BC_048")
                           ) %>% factor(., levels = c("pre","post"))
BCG$coarse_celltype_w_BCG = paste0(BCG$coarse_celltype,"_",BCG$pre_pos)
BCG <- NormalizeData(BCG)

celltypes = c("B","Treg","CD8_T","CD4_Tconv","Mono_Mac","DC","NK_ILC")

for(celltype in celltypes){
  Idents(BCG) = "coarse_celltype"
  
  print(paste0(celltype, "_start!"))
  Idents(BCG) = "coarse_celltype_w_BCG"
  deg <- Seurat::FindMarkers(object = BCG,
                             only.pos = F, 
                             logfc.threshold = 0,
                             min.pct = 0,
                             ident.1 = paste0(celltype,"_post"),
                             ident.2 = paste0(celltype,"_pre"),
                             slot = "data",
                             test.use = "MAST"
                             )
  write.table_n_2(deg, "Gene", paste0("Reproducibility/Results/Differential/BCG/UC_DOGMA_DEG_MAST_BCG_", celltype,".txt"))

  print(paste0(celltype, "_finished!"))
}


