####################################################################################
## Load the BayesPrism package
####################################################################################

source("~/my.source.R")
options(stringsAsFactors=F)
suppressWarnings(library(BayesPrism))
suppressWarnings(library(Seurat))
suppressWarnings(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))
library(ggpubr)
library(corrplot)
set.seed(1234)
setwd(path_to_wd)

####################################################################################
## Load the dataset
####################################################################################
## Extract raw count data sample*gene
bk.dat = fread_n("Reproducibility/Data/External_cohort/TCGA_BLCA_RNAseq_raw_counts.txt") %>% as.matrix() %>% t() 

scRNA = readRDS("Reproducibility/Data/Seurat/UC_DOGMA_seurat_obj_Global_RNA_HQ.rds")

# Remove post-BCG samples and Normal urothelium
Idents(scRNA) = "sample"
scRNA = subset(scRNA, idents = c('BC_027','BC_032',"BC_039","BC_040",'BC_043',"BC_044","BC_047","BC_048"), invert = TRUE)

Idents(scRNA) = 'celltype'
scRNA = subset(scRNA, idents = c('Normal'), invert = TRUE)

scRNA$cell_type = fct_collapse(scRNA$celltype,
                   Tumor = c("LUM","NRP","SQM","MES","NEC"),
                   Endothelial = c("Endothelial"),
                   B = c("Atypical_B","B_memory","B_naive",'GC_B'),
                   Plasma = c("Plasma"),
                   TNK = c("CD4_Tn","CD4_Tcm","CD4_Tsen","CD4_T_CD26","CD4_Th17","CD4_CTL","CD4_Tfh-like","CD4_Tph-like",
                           "CD4_T_proliferative","Treg_naive","Treg_effector","CD8_Tn","CD8_Tcm","CD8_Tem","CD8_Temra","CD8_Trm",
                           "CD8_Tex_1","CD8_Tex_2","CD8_T_proliferative",
                           "NK_CD56_CD49a_Hi_CD103_Hi","NK_CD56_CD49a_Hi_CD103_Lo","NK_CD56_CD49a_Lo","NK_CD56_dim","ILC3","MAIT"),
                   cDC1 = c("cDC1"),
                   mregDC = c('mregDC'),
                   MoMac_cDC2 = c("Mono","MDSC-like","TAM_TREM2","TAM_FOLR2","cDC2",'preDC'),
                   pDC = c('pDC'),
                   Mast = c('Mast'),
                   proCAF = c("proCAF"),
                   iCAF = c("iCAF_CD321",'iCAF_SLC14A1',"matCAF"),
                   CAP = c("matCAP","contCAP","vSMC")
                   )

scRNA$cell_type <- as.character(scRNA$cell_type)
scRNA$cell_state <- as.character(scRNA$cell_type)
scRNA$tumor_cluster = paste0('Tumor_', scRNA$scvi_cluster)

scRNA$cell_state[scRNA$cell_type %in% c('Tumor')] <- as.character(scRNA$tumor_cluster[scRNA$cell_type %in% c('Tumor')])

celltype_counts <- table(scRNA$cell_state)
valid_clusters <- names(celltype_counts[celltype_counts >= 100])
scRNA <- subset(scRNA, subset = cell_state %in% valid_clusters)

sc.dat = scRNA[["RNA"]]@counts %>% t() %>% as.matrix()

# Labels
cell.state.labels = scRNA@meta.data$cell_state
names(cell.state.labels) = rownames(scRNA@meta.data)

cell.type.labels = scRNA@meta.data$cell_type
names(cell.type.labels) = rownames(scRNA@meta.data)

sort(table(cell.type.labels))
table(cbind.data.frame(cell.state.labels, cell.type.labels))

####################################################################################
## QC
####################################################################################

plot.cor.phi(input=sc.dat,
             input.labels=cell.state.labels,
             title="cell state correlation",
             pdf.prefix="Reproducibility/Results/BayesPrism/1-1_QC_of_cell_state_labels", 
             cexRow=0.4, cexCol=0.4,
             margins=c(2,2))

plot.cor.phi(input=sc.dat, 
             input.labels=cell.type.labels, 
             title="cell type correlation",
             pdf.prefix="Reproducibility/Results/BayesPrism/1-2_QC_of_cell_type_labels",
             cexRow=0.5, cexCol=0.5,
             )

####################################################################################
## Filter outlier genes
####################################################################################

sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE, #return the data used for plotting. 
  pdf.prefix="Reproducibility/Results/BayesPrism/2-1_QC_filter_genes"  # "_scRNA_outlier.pdf"
)

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
    sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE,
  pdf.prefix="Reproducibility/Results/BayesPrism/2-2_QC_filter_genes"  # "_bulk_outlier.pdf"
)

# Filter outlier genes from scRNA-seq data
# Note that when sex is not identical between the reference and mixture, we recommend excluding genes from chrX and chrY. 
# We also remove lowly transcribed genes, as the measurement of transcription of these genes tend to be noise-prone. 

sc.dat.filtered <- cleanup.genes(input=sc.dat,
                                 input.type="count.matrix",
                                 species="hs", 
                                 gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                 exp.cells=5)

# Note this function only works for human data. For other species, you are advised to make plots by yourself.
plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                bulk.input = bk.dat,
                pdf.prefix="Reproducibility/Results/BayesPrism/2-3_bluk_vs_sc")

sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered,
                                       gene.type = "protein_coding")


# Select marker genes (Optional)
# performing pair-wise t test for cell states from different cell types
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],              # filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1,                                 # a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50,                             # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=32                                        # number of threads
                              )

# To subset our count matrix over the signature genes
sc.dat.filtered.pc.sig <- select.marker(sc.dat=sc.dat.filtered.pc,
                                        stat=diff.exp.stat,
                                        pval.max=0.01,
                                        lfc.min=0.1)

####################################################################################
## Construct a prism object
####################################################################################
# When using scRNA-seq count matrix as the input (recommended), user needs to specify input.type = "count.matrix". 
# The other option for input.type is "GEP" (gene expression profile) which is a cell state by gene matrix. This option is used when using reference derived from other assays, such as sorted bulk data.
# The parameter key is a character in cell.type.labels that corresponds to the malignant cell type. 
# Set to NULL if there are no malignant cells or the malignant cells between reference and mixture are from matched samples, in which case all cell types will be treated equally.

myPrism <- new.prism(
  reference=sc.dat.filtered.pc.sig, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key="Tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

####################################################################################
## Run BayesPrism
####################################################################################

bp.res <- run.prism(prism = myPrism, n.cores=32)

####################################################################################
## Extract result
####################################################################################
# extract posterior mean of cell type fraction theta
theta <- get.fraction(bp=bp.res,
                      which.theta="final",
                      state.or.type="type")

write.table_n_2(theta %>% as.data.frame(), "ID", "Reproducibility/Results/BayesPrism/bayesprism_TCGA_prop_df.txt")

save(bp.res, file="Reproducibility/Results/BayesPrism/bayesprism_result.rds")