#####################################################
# DOGMA analysis - CD4_T -
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
suppressMessages(library(CellChat))

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
suppressMessages(library(SCpubr))

suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsci))
suppressMessages(library(colorspace))
suppressMessages(library(BuenColors))

suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(escape))
suppressMessages(library(forestmodel))

suppressMessages(library(dplyr))
suppressMessages(library(jsonlite))

set.seed(1234)

#******************************************
data_dir = "Reproducibility/Data"
lineage = 'CD4_T'

DOGMA = file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")) %>% readRDS()
celltype_list = c("CD4_Tn","CD4_Tcm",'CD4_Tsen',"CD4_T_CD26","CD4_CTL","CD4_Th17","CD4_Tfh-like",
                  "CD4_Tph-like","CD4_T_proliferative","Treg_naive","Treg_effector")

#################################################
# Fig.S7C  Signature score heatmap (Primary)
#################################################
sig_df = fread_n('Reproducibility/Results/VISIONR/UC_DOGMA_CD4_T_signature_score_literature.txt')
DOGMA[["signature"]]   = CreateAssayObject(data = sig_df %>% t() %>% as.matrix())

Idents(DOGMA) = 'STAGE'
DOGMA_primary = subset(DOGMA, ident = c("post_BCG"), invert = TRUE)

Idents(DOGMA_primary) = "celltype"
DOGMA_primary = subset(DOGMA_primary, ident = c("CD4_T_proliferative"), invert = TRUE)
pseudobulk = AverageExpression(DOGMA_primary, assay = "signature", slot = "data", group.by = "celltype", return.seurat = TRUE)
pseudobulk = ScaleData(pseudobulk, scale.max = 10)  ## Default = 10

pal <- diverging_hcl(5, palette = "Cyan-magenta")
col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), pal)

average_df = pseudobulk[['signature']]$scale.data
row_order = c("Naive","Cytotoxicity","Treg","Exhaustion","Neoantigen-reactive",
              "TCR-signaling","JAK-STAT-signaling","MAPK-signaling","NFKB-signaling","IFN-signaling",
              "Cytokine","Chemokine","Costimulatory-molecules","Oxidative-phosphorylation","Glycolysis",
              "Lipid-metabolism")

ht = Heatmap(
  average_df[row_order,],
  name = "Expression",
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = TRUE, 
  show_column_names = TRUE,
  col = col_fun,
  row_names_side = "left", 
  column_names_side = "top",
  heatmap_legend_param = list(
    title = "Expression",
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2")
  )
)

paste0('Reproducibility/Results/Plots/CD4_T/FigureS7C.pdf') %>% pdf(.,w = 5, h = 5)
 draw(ht)
dev.off()


#################################################
# Fig.4C  Milo beeswarm plot
#################################################
res1 = fread_n("Reproducibility/Results/Milo/output/Milo_CD4_T_design_Organ_STAGE_contrasts_Early_H_vs_Early_L_result.txt")
res2 = fread_n("Reproducibility/Results/Milo/output/Milo_CD4_T_design_Organ_STAGE_contrasts_Advanced_vs_Early_L_result.txt")

p1 = res1 %>%
     dplyr::filter(nhood_annotation_frac > 0.5) %>%
     mutate(signif=ifelse(SpatialFDR < 0.1, logFC, 0)) %>%
     group_by(nhood_annotation)%>%
     mutate(mean_lfc_val = ifelse(!is.na(signif), logFC, 0)) %>%
     mutate(mean_lfc = mean(mean_lfc_val)) %>%
     ungroup() %>%
     dplyr::arrange(mean_lfc) %>%
     mutate(nhood_annotation=factor(nhood_annotation, 
            levels=rev(celltype_list))) %>%
     ggplot(aes(nhood_annotation, logFC)) + 
     # Grey dots first (signif = 0)
     ggbeeswarm::geom_quasirandom(data = . %>% dplyr::filter(signif == 0), 
                                  color = "grey", size = 2) +
     # Colored dots on top (signif ≠ 0)
     ggbeeswarm::geom_quasirandom(data = . %>% dplyr::filter(signif != 0), 
                                  aes(color = signif), size = 2) +
     coord_flip() +
     scale_color_gradient2(high = '#8E063B', mid = 'grey', low = '#023FA5', 
                           name = 'DA logFC\n(10% SpatialFDR)',
                           limits = c(-8, 8)) +
     theme_bw(base_size = 18) +
     geom_hline(yintercept = 0, linetype = 2) +
     xlab('') + 
     ylim(-8, 8) +
     labs(title='NMI_H_vs_NMI_L')

p2 = res2 %>%
     dplyr::filter(nhood_annotation_frac > 0.5) %>%
     mutate(signif=ifelse(SpatialFDR < 0.1, logFC, 0)) %>%
     group_by(nhood_annotation)%>%
     mutate(mean_lfc_val = ifelse(!is.na(signif), logFC, 0)) %>%
     mutate(mean_lfc = mean(mean_lfc_val)) %>%
     ungroup() %>%
     dplyr::arrange(mean_lfc) %>%
     mutate(nhood_annotation=factor(nhood_annotation, 
            levels=rev(celltype_list))) %>%
     ggplot(aes(nhood_annotation, logFC)) + 
     # Grey dots first (signif = 0)
     ggbeeswarm::geom_quasirandom(data = . %>% dplyr::filter(signif == 0), 
                                  color = "grey", size = 2) +
     # Colored dots on top (signif ≠ 0)
     ggbeeswarm::geom_quasirandom(data = . %>% dplyr::filter(signif != 0), 
                                  aes(color = signif), size = 2) +
     coord_flip() +
     scale_color_gradient2(high = '#8E063B', mid = 'grey', low = '#023FA5', 
                           name = 'DA logFC\n(10% SpatialFDR)',
                           limits = c(-8, 8)) +
     theme_bw(base_size = 18) +
     geom_hline(yintercept = 0, linetype = 2) +
     xlab('') + 
     ylim(-8, 8) +
     labs(title='MI_vs_NMI_L')

paste0("Reproducibility/Results/Plots/CD4_T/Figure4C.pdf") %>% pdf(.,width = 14, height = 5)
 p1+p2
dev.off()


#################################################
# Fig.4D  Treg fate probability
#################################################

df <- read.csv("Reproducibility/Results/CellRank2/CD4_T/fate_probabilities_dpt_Treg.csv",
               row.names = 1, check.names = FALSE)
prob_score  <- df[, c("Treg_effector_2", "Treg_effector_3")]
prob_score2 <- prob_score - 0.5   # center at 0

DOGMA$Treg_effector_2_centered <- NA_real_
overlap <- intersect(rownames(df), colnames(DOGMA))
DOGMA$Treg_effector_2_centered[overlap] <- prob_score2[overlap, "Treg_effector_2"]

p <- SCpubr::do_FeaturePlot(
  sample              = DOGMA,
  features            = "Treg_effector_2_centered",
  enforce_symmetry    = TRUE,
  diverging.palette   = "RdYlGn",
  diverging.direction = -1,
  na.value            = "grey95",   # ← brighter background
  use_viridis         = FALSE,
  plot_cell_borders   = TRUE,
  pt.size             = 4,
  border.size         = 1,
  raster              = TRUE,
  raster.dpi          = 1024
)

paste0("Reproducibility/Results/Plots/CD4_T/Figure4D_fate_prob_UMAP_w_bg.pdf") %>% pdf(w=3, h=5)
 plot(p)
dev.off()

#################################################
# Fig.S7E  eRegulon score box plot
#################################################
eRegulon_df = fread_n("Reproducibility/Results/scenicplus/CD4_T/eRegulon_df.txt")
eRegulon_meta = fread_n("Reproducibility/Results/scenicplus/CD4_T/eRegulon_metadata.txt")
colnames(eRegulon_meta) = c('celltype','celltype_ATAC')

TFs = c("MAF_extended_+/+_(106g)","RUNX1_extended_+/+_(638g)","TBX21_extended_+/+_(29g)")
celltype_keep <- c("CD4_Tn","CD4_Tcm",'CD4_Tsen',"CD4_T_CD26","CD4_CTL","CD4_Th17","CD4_Tfh-like","CD4_Tph-like")
celltype_colors <- c(
  "CD4_Tn"       = '#1f77b4',
  "CD4_Tcm"      = '#ff7f0e',
  "CD4_Tsen"     = '#279e68',
  "CD4_T_CD26"   = '#d62728',
  "CD4_CTL"      = '#aa40fc',
  "CD4_Th17"     = '#8c564b',
  "CD4_Tfh-like" = '#e377c2',
  "CD4_Tph-like" = '#b5bd61'
)

plot_df <- eRegulon_df %>%
  mutate(celltype = eRegulon_meta$celltype) %>%
  dplyr::select(all_of(TFs), celltype) %>%
  dplyr::filter(celltype %in% celltype_keep) %>%
  pivot_longer(cols = all_of(TFs), names_to = "TF", values_to = "score") %>%
  mutate(
    celltype = factor(celltype, levels = celltype_keep),
    TF = factor(TF, levels = TFs)
  )

# Headroom for significance bars
y_max <- max(plot_df$score, na.rm = TRUE)
y_expand_top <- 0.15    # 15% headroom to avoid crowding of significance bars

p <- ggplot(plot_df, aes(x = celltype, y = score)) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  geom_quasirandom(
    aes(color = celltype),
    method = "pseudorandom",
    size = 2.2, alpha = 0.45, shape = 16
  ) +
  facet_wrap(~ TF, ncol = 1, scales = "free_y") +   # <- independent Y scales
  scale_color_manual(values = celltype_colors) +
  labs(x = "Cell Type", y = "eRegulon score") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

pdf("Reproducibility/Results/Plots/CD4_T/FigureS7E.pdf", width = 7, height = 9)
 print(p)
dev.off()


#########################################################
## Fig.S7D/F  Gene track with co-accessible regions
#########################################################
source("Reproducibility/Scripts/Source/Seurat_source.R")

DOGMA$celltype = factor(DOGMA$celltype, levels=celltype_list)

# CD4_Tconv
region_list <- list(
  IL17A = c('chr6-52165000-52230000')
)

motif_list = list(
  IL17A = c('RORC','RORA')
)

cols <- c('#1f77b4','#ff7f0e','#279e68','#d62728','#aa40fc','#8c564b','#e377c2','#b5bd61')

Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='CD4_T',group_list=celltype_keep,gene='IL17A',
               region_list=region_list, motif_list=motif_list,cutoff=0.3,cols=cols,
               path='Reproducibility/Results/Plots/CD4_T/FIgureS7D.pdf')

# Treg
celltype_list=c("Treg_naive","Treg_effector")
region_list <- list(
  CCR8 = c('chr3-39308000-39342000')
)

motif_list = list(
  CCR8 = c('BATF','SOX4','REL')
)

cols <- c('#aec7e8', '#ffbb78')

Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='CD4_T',group_list=celltype_list,gene='CCR8',
               region_list=region_list, motif_list=motif_list,cutoff=0.3,cols=cols,
               path='Reproducibility/Results/Plots/CD4_T/FIgureS7F.pdf')



#########################################################
## Fig.4G  Gene track with co-accessible regions
#########################################################

df <- read.csv("Reproducibility/Results/CellRank2/CD4_T/fate_probabilities_dpt_Treg.csv", row.names = 1)
prob_score = df[, c("Treg_effector_2", "Treg_effector_3")] 

# Get barcodes
top_barcodes <- dplyr::arrange(prob_score, desc(Treg_effector_2)) %>% head(2000) %>% rownames()
bottom_barcodes <- dplyr::arrange(prob_score, desc(Treg_effector_3)) %>% head(2000) %>% rownames()

###################################
Idents(DOGMA) = "celltype"
Treg = subset(DOGMA, ident= c('Treg_naive', 'Treg_effector'))
Treg$state <- ifelse(Treg$celltype == "Treg_naive", "naive",
              ifelse(colnames(Treg) %in% top_barcodes, "State1",
              ifelse(colnames(Treg) %in% bottom_barcodes, "State2", "Other"))) %>% 
              factor(, levels = c('naive','State1', 'State2','Other'))
group_list=c('naive','State1', 'State2')

region_list <- list(
  CXCR3 = c('chrX-71610000-71630000')
)

motif_list = list(
  CXCR3 = c('TBX21')
)

cols <- c("#756bb1","#F46D43","#66BD63")

DefaultAssay(Treg)='ATAC'
Plot_genetrack(Obj=Treg,group.by='state',lineage='CD4_T',group_list=group_list,gene='CXCR3',
               region_list=region_list, motif_list=motif_list,cutoff=0.1,cols=cols,
               path='Reproducibility/Results/Plots/CD4_T/FIgure4G.pdf')


#########################################################
## Fig.S7H  Heatmap of functional genes
#########################################################

treg_nTreg_genes <- c("TCF7","BCL6")
treg_core_signature <- c("FOXP3","IL2RA","IL2RB","IKZF2","CTLA4","GPR83",'IL10')
eTreg_tissue_adaptation <- c("TRAF1","TRAF4","NFKB1","NFKB2","RELB","NFKBIZ","CREM","IL1R2",
                             "TNFRSF4","TNFRSF9","TNFRSF18","CD83",'CCR8')
Antigen_processing = c("HLA-DRA","HLA-DRB1")
isg_core <- c(
  "STAT1","STAT2","IRF7",
  "ISG15","MX1",
  "IFIT1","IFI6","IFI44L","CXCL10")
genes_stress <- c("HSPA1A","HSPA1B","HSP90AA1","HMOX1","GADD45B","TXN","SESN2","DNAJB1")

pseudobulk_tmp = AverageExpression(Treg, assay = "RNA", slot = "data", group.by = "state", return.seurat = TRUE)
pseudobulk_tmp = ScaleData(pseudobulk_tmp, scale.max = 10)
RNA_scaled_df = pseudobulk_tmp[['RNA']]$scale.data %>% 
                as.data.frame() %>% dplyr::select(., group_list)

palette = rev(c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac"))
col_fun <- colorRamp2(c(-2,-4/3,-2/3, 0, 2/3,4/3, 2), palette)

# Create the heatmap
ht1 <- Get_Heatmap(df = RNA_scaled_df, genes_list = treg_nTreg_genes, colors = col_fun)
ht2 <- Get_Heatmap(df = RNA_scaled_df, genes_list = treg_core_signature, colors = col_fun)
ht3 <- Get_Heatmap(df = RNA_scaled_df, genes_list = eTreg_tissue_adaptation, colors = col_fun)
ht4 <- Get_Heatmap(df = RNA_scaled_df, genes_list = Antigen_processing, colors = col_fun)
ht5 <- Get_Heatmap(df = RNA_scaled_df, genes_list = isg_core, colors = col_fun)
ht6 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_stress, colors = col_fun)

paste0('Reproducibility/Results/Plots/CD4_T/FigureS7H.pdf') %>% pdf(.,width = 13, height = 2)
 draw(ht1+ht2+ht5+ht6+ht3+ht4)
dev.off()


#########################################################
## Fig.S7I  Differential analysis between Treg states
#########################################################
# RNA
data = fread_n("Reproducibility/Results/Differential/UC_DOGMA_DEG_MAST_Treg_by_pseudotime_state.txt")
data$avg_log2FC[data$avg_log2FC > 4] <- 4
data$avg_log2FC[data$avg_log2FC < -4] <- -4
label_up = c("ISG15", "IFI44L", "MX1", "CXCR3", 'FOXP3','TBX21')
label_down = c("TNFRSF4", "TNFRSF8", "TNFRSF9","RELB",'NFKB1','IL2RA')

paste0("Reproducibility/Results/Plots/CD4_T/FigureS7I.pdf") %>% pdf(.,width=3,height=3)
  Volcano_DEG_label_v2(data = data, celltype = "state", label_up = label_up, label_down = label_down, log2FC = 0.5, ylim = 150, thresh = 1e-150)
dev.off()


#########################################################
## Fig.S7K  Scatterplot of neighbourhoods
#########################################################

df = fread_n('Reproducibility/Results/Milo/output/Milo_CD4_T_design_Organ_STAGE2_contrasts_MI_vs_NMI_result.txt')
df$nhood_broad_annotation = fct_collapse(df$nhood_annotation,
                CD4_Tconv = c("CD4_Tn","CD4_Tcm",'CD4_Tsen',"CD4_T_CD26","CD4_CTL","CD4_Th17",
                              "CD4_Tfh-like","CD4_Tph-like","CD4_T_proliferative"),
                Treg = c("Treg_effector","Treg_naive"))

p1 = df %>%
    dplyr::filter(nhood_broad_annotation == 'Treg') %>%
    ggplot(aes(logFC, nhood_IFN_signaling)) +
    geom_point(size=1, alpha=0.5, shape=16) +
    facet_wrap(nhood_broad_annotation~.) +
    guides(fill=guide_legend(title='', override.aes = list(size=2), ncol=1)) +
    theme_bw(base_size=24) +
    scale_fill_brewer(palette='Spectral', name='') +
    geom_vline(xintercept=0, linetype=2) +
    xlab("DA logFC") + ylab("IFN signature") +
    theme(legend.position='top') +
    ggpubr::stat_cor(size=6) +
    xlim(-7, 7) + ylim(0,1.4) 

p2 = df %>%
    dplyr::filter(nhood_broad_annotation == 'Treg') %>%
    ggplot(aes(logFC, nhood_Glycolysis)) +
    geom_point(size=1, alpha=0.5, shape=16) +
    facet_wrap(nhood_broad_annotation~.) +
    guides(fill=guide_legend(title='', override.aes = list(size=2), ncol=1)) +
    theme_bw(base_size=24) +
    scale_fill_brewer(palette='Spectral', name='') +
    geom_vline(xintercept=0, linetype=2) +
    xlab("DA logFC") + ylab("IFN signature") +
    theme(legend.position='top') +
    ggpubr::stat_cor(size=6) +
    xlim(-7, 7) + ylim(0,2)

pdf("Reproducibility/Results/Plots/CD4_T/FigureS7K.pdf", width = 4, height = 10)
 print(p1/p2)
dev.off()


#########################################################
## Fig.4G  Differential analysis in eTregs
#########################################################
# RNA
data = fread_n("Reproducibility/Results/Differential/UC_DOGMA_DEG_MAST_Treg_effector.txt")
data$avg_log2FC[data$avg_log2FC > 4] <- 4
data$avg_log2FC[data$avg_log2FC < -4] <- -4
label_up = c('IL12RB2','IFI6','ICOS','CXCR6','IFI44L','TBX21','CXCR3','FOXP3')
label_down = c('REL','BCL2','NFKB1','SEMA4A','NR4A2')

paste0("Reproducibility/Results/Plots/CD4_T/Figure4G_RNA.pdf") %>% pdf(.,width=3,height=3)
  Volcano_DEG_label_v2(data = data, celltype = "MI enr nhoods", label_up = label_up, label_down = label_down, log2FC = 0.5)
dev.off()

###################################
# ADT
data2 = fread_n("Reproducibility/Results/Differential/UC_DOGMA_DEG_MAST_Treg_effector_ADT.txt")
data2$avg_log2FC[data2$avg_log2FC > 0.7] <- 0.7
data2$avg_log2FC[data2$avg_log2FC < -0.7] <- -0.7

rownames(data2) = take_factor(rownames(data2),3,'-')
label_up = c('ICOS','CD73','PD1')
label_down = c('CD71',"CD177","CD27")

paste0("Reproducibility/Results/Plots/CD4_T/Figure4G_ADT.pdf") %>% pdf(.,width=3,height=3)
  Volcano_DEG_label_v2(data = data2, celltype = "MI enr nhoods", label_up = label_up, label_down = label_down, ylim = 60, 
                        log2FC = 0.25, thresh = 1e-50)
dev.off()

###################################
# TF activity
data3 = fread_n("Reproducibility/Results/Differential/UC_DOGMA_DEG_wilcox_Treg_effector_TF.txt")
label_up = c('FOXP3','TBX21','CTCF','ETS1','IRF4','STAT1','IRF3')
label_down = c('RUNX1','JUN', 'STAT4')
paste0("Reproducibility/Results/Plots/CD4_T/Figure4G_TF.pdf") %>% pdf(.,width=3,height=3)
  Volcano_DEG_label_v2(data = data3, celltype = "MI enr nhoods", label_up = label_up, label_down = label_down, ylim = 100, log2FC = 0.5)
dev.off()

#########################################################
## Fig.4H/4I/S7L  CCI analysis (CellChat)
#########################################################
scRNA = paste0("Reproducibility/Data/Seurat/UC_DOGMA_seurat_obj_Global_RNA_HQ.rds") %>% readRDS()

celltype_list = c("CD4_Tn","CD4_Tcm",'CD4_Tsen',"CD4_T_CD26","CD4_CTL","CD4_Th17","CD4_Tfh-like",
                  "CD4_Tph-like","CD4_T_proliferative","Treg_naive","Treg_effector")

scRNA$coarse_celltype2 = fct_collapse(scRNA$celltype, 
    MSC = c('proCAF','iCAF_SLC14A1','iCAF_CD321','matCAF','matCAP','contCAP','vSMC'),
    Endothelial = c('Endothelial'),
    B = c('Atypical_B','B_memory','B_naive','GC_B'),
    Plasma = c('Plasma'),
    CD4_Tconv = c("CD4_Tn","CD4_Tcm",'CD4_Tsen',"CD4_T_CD26","CD4_CTL","CD4_Th17",
                  "CD4_Tfh-like","CD4_Tph-like","CD4_T_proliferative"),
    CD8_T = c('CD8_T_proliferative','CD8_Tn','CD8_Tcm','CD8_Tem','CD8_Temra',
              'CD8_Trm','CD8_Tex_1','CD8_Tex_2'),
    NK_ILC = c('ILC3','MAIT','NK_CD56_CD49a_Hi_CD103_Hi','NK_CD56_CD49a_Hi_CD103_Lo','NK_CD56_CD49a_Lo','NK_CD56_dim'),
    Mono = c('MDSC-like','Mono'),
    TAM_cDC2 = c('TAM_FOLR2','TAM_TREM2','cDC2')
)
scRNA$samples = scRNA$sample
scRNA$samples = factor(scRNA$samples)

Idents(scRNA) = "coarse_celltype2"
scRNA_tmp = subset(scRNA, idents = c('B','Plasma','MSC','CD4_Tconv','Treg_naive','Treg_effector',
          'CD8_T','NK_ILC','Mono','TAM_cDC2','cDC1','mregDC','pDC'))
scRNA_tmp$coarse_celltype2 = scRNA_tmp$coarse_celltype2[,drop=TRUE]

###############################################################################
Idents(scRNA_tmp) = 'STAGE'

for(tmp_STAGE in c('Early_L','Early_H','Advanced')){
  print(paste0(tmp_STAGE, ' start'))

  scRNA_tmp_STAGE = subset(scRNA_tmp, idents = tmp_STAGE)

  # 1) Data input & processing and initialization of CellChat object
  cellchat <- createCellChat(object = scRNA_tmp_STAGE, group.by = "coarse_celltype2", assay = "RNA")
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB)
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # 2) Inference of cell-cell communication network
  ptm = Sys.time()
  cellchat@idents <- droplevels(cellchat@idents)

  options(future.globals.maxSize = 16 * 1024^3)
  future::plan("multisession", workers = 4)
  print(future::nbrOfWorkers())
  print(getOption("future.globals.maxSize")) 

  cellchat <- computeCommunProb(cellchat, type = "truncatedMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  # 3) Extract the inferred cellular communication network as a data frame
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  saveRDS(cellchat, file = paste0("Reproducibility/Results/CellChat/UC_DOGMA_STAGE_",tmp_STAGE,".rds"))

  print(paste0(tmp_STAGE, ' finish'))
}

cellchat.H <- readRDS("Reproducibility/Results/CellChat/UC_DOGMA_STAGE_Early_H.rds")
cellchat.A <- readRDS("Reproducibility/Results/CellChat/UC_DOGMA_STAGE_Advanced.rds")

cellchat.H <- netAnalysis_computeCentrality(cellchat.H, slot.name = "netP")
cellchat.A <- netAnalysis_computeCentrality(cellchat.A, slot.name = "netP")

object.list <- list(H = cellchat.H, A = cellchat.A)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

###############################################################################
## Visually compare cell-cell communication using Circle plot

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
paste0('Reproducibility/Results/Plots/CD4_T/FIgure4H.pdf') %>% pdf_3(., w=12, h=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pathways.show <- c("IFN-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
paste0('Reproducibility/Results/Plots/CD4_T/FIgureS7L.pdf') %>% pdf_3(., w=12, h=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

###############################################################################
## Identify the up/down-regulated signaling ligand-receptor pairs

pos.dataset = "A"
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "A",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "H",ligand.logFC = -0.05, receptor.logFC = NULL)

pairLR.use.up = dplyr::filter(net.up, pathway_name %in% c('CCL','CXCL'))
pairLR.use.up = pairLR.use.up[, "interaction_name", drop = F] 

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(5,4,9), targets.use = 12, 
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

paste0('Reproducibility/Results/Plots/CD4_T/FIgure4I.pdf') %>% pdf_3(., w=4, h=3)
 gg1
dev.off()