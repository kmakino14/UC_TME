#####################################################
# DOGMA analysis - Myeloid -
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
suppressMessages(library(escape))
suppressMessages(library(forestmodel))

suppressMessages(library(dplyr))
suppressMessages(library(jsonlite))

set.seed(1234)

#******************************************
data_dir = "Reproducibility/Data"
lineage = 'Myeloid'

DOGMA = file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")) %>% readRDS()

#################################################
# Fig.S7C  Functional genes heatmap
#################################################
celltype = c('Mono','MDSC-like','TAM_TREM2','TAM_FOLR2')

Idents(DOGMA) = 'celltype'
DOGMA_part = subset(DOGMA, ident = celltype)

pseudobulk_tmp = AverageExpression(DOGMA_part, assay = "RNA", slot = "data", group.by = "celltype", return.seurat = TRUE)
pseudobulk_tmp$celltype = colnames(pseudobulk_tmp) %>% 
                          factor(., levels = celltype)
pseudobulk_tmp = ScaleData(pseudobulk_tmp, scale.max = 10)

RNA_scaled_df = pseudobulk_tmp[['RNA']]@scale.data

p = rev(c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac"))
col_fun <- colorRamp2(c(-2,-4/3,-2/3, 0, 2/3,4/3, 2), p)

ht <- function(df, max, min, col ,cluster_columns = TRUE){
    ht = Heatmap(
         df,                              # Transpose to match row/column orientation
         name = "Expression",             # Name of the heatmap
         cluster_rows = FALSE,            # Enable row clustering (equivalent to dendrogram=True)
         cluster_columns = cluster_columns,       # Enable column clustering (equivalent to dendrogram=True)
         show_row_names = TRUE,           # Display row names (signatures/genes)
         show_column_names = TRUE,        # Display column names (cell types)
         col = col,                       # Apply the color function
         row_names_side = "left",         # Position row names on the left
         column_names_side = "top",       # Position column names on the top
         heatmap_legend_param = list(
            title = "Expression",
            at = c(min, 0, max),              # Set legend ticks similar to vmin and vmax
            labels = c(as.character(min), "0", as.character(max))
         )
         )
    return(ht)
}

# Create the heatmap
ht1 <- ht(t(RNA_scaled_df[c('S100A6','TREM1','FYN'),celltype]),max=2,min=-2,col=col_fun) 
ht2 <- ht(t(RNA_scaled_df[c('CEBPB','VEGFA','NCF2'),celltype]),max=2,min=-2,col=col_fun) 
ht3 <- ht(t(RNA_scaled_df[c('APOC1','C1QA','APOE','CD163','MRC1'),celltype]),max=2,min=-2,col=col_fun) 
ht4 <- ht(t(RNA_scaled_df[c('TNF','CD80','CD86','CD40','CCL5'),celltype]),max=2,min=-2,col=col_fun)
ht5 <- ht(t(RNA_scaled_df[c('FCGR2B','FN1','MARCO','SPP1','ARG2','CSF1R','MSR1','CCL4','CD276'),celltype]),max=2,min=-2,col=col_fun)
ht6 <- ht(t(RNA_scaled_df[c('IL6','CXCL2','CCL3','CCL4','CXCL8','CXCL1','CCL2','CXCL3'),celltype]),max=2,min=-2,col=col_fun) 
ht7 <- ht(t(RNA_scaled_df[c('LIPA','C1QA','C1QB','C1QC','GPNMB','TREM2','MERTK'),celltype]),max=2,min=-2,col=col_fun) 
ht8 <- ht(t(RNA_scaled_df[c('S100A9','S100A8','VCAN','AREG','THBS1','S100A12','FCN1','CD44'),celltype]),max=2,min=-2,col=col_fun) 
ht9 <- ht(t(RNA_scaled_df[c('FOLR2','PLTP','LYVE1'),celltype]),max=2,min=-2,col=col_fun)
ht10 <- ht(t(RNA_scaled_df[c('RGS2','HSP90AA1','DNAJB1','HSP90AB1','ZFAND2A','BAG3','HSPB1','HSPA1A','SOCS3'),celltype]),max=2,min=-2,col=col_fun)

paste0('Reproducibility/Results/Plots/Myeloid/FigureS7C.pdf') %>% pdf(.,w = 15, h = 3)
 draw(ht1+ht2+ht3+ht4+ht5+ht6+ht7+ht8+ht9+ht10)
dev.off()


#################################################
# Fig.S7D  Signature violin plot
#################################################

sig_df = fread_n('Reproducibility/Results/VISIONR/UC_DOGMA_TAM_signature_score_literature.txt') 
DOGMA_part[["signature"]]   = CreateAssayObject(data = sig_df %>% t() %>% as.matrix())
DOGMA_part$celltype = factor(DOGMA_part$celltype, levels = celltype)
cols = c("#2372A9","#EF7C21",'#2D9865','#CA2A28')

vlnplot_w_box(Seurat = DOGMA_part,
              assay = 'signature', 
              slot = 'data', 
              features = rownames(DOGMA_part[["signature"]]), 
              file_name = "Reproducibility/Results/Plots/Myeloid/FigureS7D.pdf", 
              cols = cols)


#################################################
# Fig.S7E  Milo beeswarm plot
#################################################
res1 = fread_n("Reproducibility/Results/Milo/output/Milo_Myeloid_design_Organ_STAGE_contrasts_Early_H_vs_Early_L_result.txt")
res2 = fread_n("Reproducibility/Results/Milo/output/Milo_Myeloid_design_Organ_STAGE_contrasts_Advanced_vs_Early_L_result.txt")

p1 = res1 %>%
     dplyr::filter(nhood_annotation_frac > 0.5) %>%
     mutate(signif=ifelse(SpatialFDR < 0.1, logFC, 0)) %>%
     group_by(nhood_annotation)%>%
     mutate(mean_lfc_val = ifelse(!is.na(signif), logFC, 0)) %>%
     mutate(mean_lfc = mean(mean_lfc_val)) %>%
     ungroup() %>%
     dplyr::arrange(mean_lfc) %>%
     mutate(nhood_annotation=factor(nhood_annotation, 
            levels=c('Mast',"pDC",'preDC','mregDC','cDC2','cDC1','TAM_FOLR2','TAM_TREM2','MDSC-like','Mono'))) %>%
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
            levels=c('Mast',"pDC",'preDC','mregDC','cDC2','cDC1','TAM_FOLR2','TAM_TREM2','MDSC-like','Mono'))) %>%
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

paste0("Reproducibility/Results/Plots/Myeloid/FigureS7E.pdf") %>% pdf(.,width = 14, height = 5)
 p1+p2
dev.off()


#################################################
# Fig.S7F  eRegulon score box plot
#################################################
eRegulon_df = fread_n("Reproducibility/Results/scenicplus/Myeloid/eRegulon_df.txt")
eRegulon_meta = fread_n("Reproducibility/Results/scenicplus/Myeloid/eRegulon_metadata.txt")
colnames(eRegulon_meta) = c('celltype','celltype_ATAC')

TFs = c("RELB_extended_+/+_(392g)", "IKZF1_direct_+/+_(388g)")
celltype_keep <- c('cDC1','cDC2','mregDC','preDC','pDC')
celltype_colors <- c(
  "cDC1"   = '#aa40fc',
  "cDC2"   = '#8c564b',
  "mregDC" = '#e377c2',
  "preDC"  = '#b5bd61',
  "pDC"    = '#17becf'
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

comparisons <- combn(levels(plot_df$celltype), 2, simplify = FALSE)

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
  facet_wrap(~ TF, ncol = 2, scales = "free_y") +   # <- independent Y scales
  scale_color_manual(values = celltype_colors) +
  labs(x = "Cell Type", y = "eRegulon score") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none") +
  ggpubr::stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "BH",
    hide.ns = TRUE,
    tip.length = 0,
    step.increase = 0.05
  )

pdf("Reproducibility/Results/Plots/Myeloid/FigureS7F.pdf", width = 10, height = 6)
 print(p)
dev.off()


#########################################################
## Fig.S7G  Gene track with co-accessible regions
#########################################################
source("Reproducibility/Scripts/Source/Seurat_source.R")

DOGMA$celltype = factor(DOGMA$celltype, levels=c("Mono","MDSC-like","TAM_TREM2","TAM_FOLR2","cDC1","cDC2","mregDC","preDC","pDC","Mast"))

# DCs
celltype_list=c('cDC1','cDC2','mregDC','preDC','pDC')
region_list <- list(
  CD86 = c('chr3-121990000-122130000')
)

motif_list = list(
  CD86 = c('REL','IKZF1')
)

cols <- c('#aa40fc',"#8c564b","#e377c2",'#b5bd61','#17becf')

Plot_genetrack(Obj=DOGMA, group.by = 'celltype', lineage='Myeloid',
               group_list=celltype_list,gene='CD86',
               region_list=region_list, motif_list=motif_list,cutoff=0.2,cols=cols,
               path='Reproducibility/Results/Plots/Myeloid/FigureS7G_CD86.pdf')

# Mono/Mac
celltype_list=c("Mono","MDSC-like","TAM_TREM2","TAM_FOLR2")
region_list <- list(
  VEGFA = c('chr6-43760000-43810000')
)

motif_list = list(
  VEGFA = c('ARNT::HIF1A','SMAD2')
)

cols <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728")

Plot_genetrack(Obj=DOGMA,group.by = 'celltype', lineage='Myeloid',
               group_list=celltype_list,gene='VEGFA',
               region_list=region_list, motif_list=motif_list,cutoff=0.2,cols=cols,
               path='Reproducibility/Results/Plots/Myeloid/FigureS7G_VEGFA.pdf')
