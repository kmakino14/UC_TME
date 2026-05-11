#####################################################
# DOGMA analysis - CD8_T_ILC -
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
suppressMessages(library(ggalluvial))
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
lineage = 'CD8_T_NK_ILC'

DOGMA = file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")) %>% readRDS()
celltype_list = c("CD8_Tn","CD8_Tcm",'CD8_Tem',"CD8_Temra","CD8_Trm","CD8_Tex_1",
        "CD8_Tex_2","CD8_T_proliferative","NK_CD56_CD49a_Hi_CD103_Hi","NK_CD56_CD49a_Hi_CD103_Lo",
        "NK_CD56_CD49a_Lo","NK_CD56_dim",'ILC3','MAIT')


#################################################
# Fig.S10D  Milo beeswarm plot
#################################################
res1 = fread_n("Reproducibility/Results/Milo/output/Milo_CD8_T_NK_ILC_design_Organ_STAGE_contrasts_Early_H_vs_Early_L_result.txt")
res2 = fread_n("Reproducibility/Results/Milo/output/Milo_CD8_T_NK_ILC_design_Organ_STAGE_contrasts_Advanced_vs_Early_L_result.txt")

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

paste0("Reproducibility/Results/Plots/CD8_T_NK_ILC/FigureS10D.pdf") %>% pdf(.,width = 18, height = 6)
 p1+p2
dev.off()


#################################################
# Fig.5J  CD8 fate probability
#################################################

df <- read.csv("Reproducibility/Results/CellRank2/CD8_T_NK_ILC/fate_probabilities_dpt_CD8.csv", row.names = 1)
prob_score = df[, c("CD8_Tex_2_2", "CD8_Temra_2")] 
prob_score2 = prob_score-0.5 

# Create Seurat object with dummy counts
seurat_obj <- CreateSeuratObject(counts = prob_score2 %>% t(), meta.data = df)

# Attach UMAP coordinates
seurat_obj[["umap"]] <- CreateDimReducObject(
  embeddings = df[, c("UMAP_1", "UMAP_2")] %>% as.matrix(),
  key = "UMAP_",
  assay = "RNA"
)

p = SCpubr::do_FeaturePlot(
  seurat_obj,
  features = 'CD8-Tex-2-2',
  enforce_symmetry = TRUE,
  diverging.palette = "PiYG",
  diverging.direction = -1,
  plot_cell_borders = TRUE,
  pt.size = 4,
  border.size = 1,
  raster = TRUE,
  raster.dpi = 1024)

paste0("Reproducibility/Results/Plots/CD8_T_NK_ILC/Figure5J_fate_prob_UMAP.pdf") %>% pdf(w=3, h=5)
 plot(p)
dev.off()


#########################################################
## Fig.S10H/S11E  Gene track with co-accessible regions
#########################################################

DOGMA$celltype = factor(DOGMA$celltype, levels=celltype_list)

# CD8_T
region_list <- list(
  CTLA4 = c('chr2-203810000-203900000'),
  LAYN = c('chr11-111528000-111545000')
)

motif_list = list(
  CTLA4 = c('PRDM1','BATF','NFATC3'),
  LAYN  = c('PRDM1','BATF','NFKB1','PPARG')
)

cols <- c('#1aafc9','#ff6f61','#6a9f58','#f2b134','#d883b7','#6b4c9a','#c0592c')
celltype_keep <- c("CD8_Tn","CD8_Tcm",'CD8_Tem',"CD8_Temra","CD8_Trm","CD8_Tex_1","CD8_Tex_2")

Plot_genetrack(Obj=DOGMA,lineage='CD8_T_NK_ILC',celltype_list=celltype_keep,gene='CTLA4',
               region_list=region_list, motif_list=motif_list,cutoff=0.4,cols=cols,
               path='Reproducibility/Results/Plots/CD8_T_NK_ILC/FIgureS10H_CTLA4.pdf')

Plot_genetrack(Obj=DOGMA,lineage='CD8_T_NK_ILC',celltype_list=celltype_keep,gene='LAYN',
               region_list=region_list, motif_list=motif_list,cutoff=0.2,cols=cols,
               path='Reproducibility/Results/Plots/CD8_T_NK_ILC/FIgureS10H_LAYN.pdf')

# NK
region_list <- list(
  KLRC1 = c('chr12-10430000-10575000'),
  ITGAE = c('chr17-3710000-3815000')
)

motif_list = list(
  KLRC1 = c('PRDM1','RUNX3'),
  ITGAE = c('ELF1','RUNX3')
)

cols <- c( '#3fb68b','#bc4b51','#7199c6','#e84d8a')
celltype_keep <- c("NK_CD56_CD49a_Hi_CD103_Hi","NK_CD56_CD49a_Hi_CD103_Lo",
                   "NK_CD56_CD49a_Lo","NK_CD56_dim")

Plot_genetrack(Obj=DOGMA,lineage='CD8_T_NK_ILC',celltype_list=celltype_keep,gene='KLRC1',
               region_list=region_list, motif_list=motif_list,cutoff=0.2,cols=cols,
               path='Reproducibility/Results/Plots/CD8_T_NK_ILC/FIggureS11E_KLRC1.pdf')

Plot_genetrack(Obj=DOGMA,lineage='CD8_T_NK_ILC',celltype_list=celltype_keep,gene='ITGAE',
               region_list=region_list, motif_list=motif_list,cutoff=0.2,cols=cols,
               path='Reproducibility/Results/Plots/CD8_T_NK_ILC/FIggureS11E_ITGAE.pdf')


#########################################################
## Fig.S11B  Functional gene heatmap
#########################################################

Idents(DOGMA) = "celltype"
DOGMA_NK = subset(DOGMA, idents = c('NK_CD56_CD49a_Hi_CD103_Hi', 'NK_CD56_CD49a_Hi_CD103_Lo', 
                                    'NK_CD56_CD49a_Lo','NK_CD56_dim'))

genes_inf = c('IFNG','TNF','IL18','IL15','IL7','IL6','IL1B','CXCL9','CXCL10','CCL5','CCL4','CCL3','CCL2')
genes_stress = c('ZFP36L1','ZFAND2A','UBC','SOCS3','SLC2A3','RGS2','JUN',
                 'HSPH1','HSPB1','HSPA6','HSPA1B','HSPA1A','HSP90AB1','HSP90AA1','DNAJB1','BAG3') 
genes_cytotoxic = c('CTSW','PRF1','GNLY','GZMK','GZMM','GZMH','GZMB','GZMA') 
genes_activate_HLA = c('CD160','KIR2DL4')
genes_activate_HLA_ind = c('HCST','FCGR3A','NCR3','NCR1','KLRK1','CRTAM')
genes_Co_receptors = c('CD226','SLAMF7','SLAMF6','CD244','TNFRSF9','CD59')
genes_inhibitory_HLA = c('KLRC1','KIR2DL1','KIR3DL2','KIR3DL1','LILRB1','LAG3')
genes_inhibitory_HLA_ind = c('KLRB1','HAVCR2','SIGLEC9','LAIR1','CD300A','TIGIT','CD96')
genes_adaptive = c('KLRC2','CD3E',"CTLA4",'ENTPD1','TOX','CD52','CCL5','IL32','PDCD1')
genes_circulatory = c('S1PR5','KLF2','CX3CR1','KLRG1')

pseudobulk_tmp = AverageExpression(DOGMA_NK, assay = "RNA", slot = "data", group.by = "celltype", return.seurat = TRUE)
pseudobulk_tmp = ScaleData(pseudobulk_tmp, scale.max = 10)
RNA_scaled_df = pseudobulk_tmp[['RNA']]$scale.data

palette = rev(c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac"))
col_fun <- colorRamp2(c(-2,-4/3,-2/3, 0, 2/3,4/3, 2), palette)

# Create the heatmap
ht1 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_inf, colors = col_fun)
ht2 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_stress, colors = col_fun)
ht3 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_cytotoxic, colors = col_fun)
ht4 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_activate_HLA, colors = col_fun)
ht5 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_activate_HLA_ind, colors = col_fun)
ht6 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_Co_receptors, colors = col_fun)
ht7 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_inhibitory_HLA, colors = col_fun)
ht8 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_inhibitory_HLA_ind, colors = col_fun)
ht9 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_adaptive, colors = col_fun)
ht10 <- Get_Heatmap(df = RNA_scaled_df, genes_list = genes_circulatory, colors = col_fun)

paste0('Reproducibility/Results/Plots/CD8_T_NK_ILC/FIggureS11B.pdf') %>% pdf(.,width = 13, height = 2)
 draw(ht1+ht2+ht3+ht4+ht5+ht6+ht7+ht8+ht9+ht10)
dev.off()


#########################################################
## Fig.5K  in silico perturbation module heatmap
#########################################################
scenic_res <- fread("Reproducibility/Results/scenicplus/CD8_T_NK_ILC/tf_to_gene_adj.tsv") %>% as.data.frame()
hs_module = fread_n("Reproducibility/Results/Hotspot/CD8_T_NK_ILC/UC_DOGMA_CD8_T_NK_ILC_Hotspot_module_df.txt") %>%
            dplyr::filter(., Module %in% c(2,5)) %>%
            dplyr::filter(., FDR<0.01)  # FDR<0.01
genes = rownames(hs_module)

# Subset by CD8_Tex genes
scenic_res2 <- scenic_res %>% dplyr::filter(target %in% genes)
scenic_res2$importance2 = scenic_res2$importance * sign(scenic_res2$rho)

ranks <- rank(abs(scenic_res2$importance2), ties.method = "average")
scenic_res2$importance2_ranknorm <- sign(scenic_res2$importance2) * (ranks / max(ranks))

df_wide <- scenic_res2 %>% 
  dplyr::select(TF, target, importance2_ranknorm) %>%  # importance
  pivot_wider(names_from = target, values_from = importance2_ranknorm, values_fill = 0) %>%
  column_to_rownames("TF")

tf_active_targets <- rowSums(df_wide > 0)

# Keep TFs with ≥10 non-zero entries
filtered_tfs <- names(tf_active_targets[tf_active_targets >= 200])
additional_filtering = c("ZNF518A", "BCLAF1", "SMARCA5", "JUND", "JUN", "CHD1", "CHD2", "HBP1", 
           "NR4A2", "CREM", "ZNF10", "STAT4", "FOSL2", "NFATC1", "ETV6", "CLOCK", 
           "MAF", "REL", "NFKB1", "NFE2L2", "TFDP2", "RORA", "BACH2")
filtered_tfs_v2 = setdiff(filtered_tfs, additional_filtering)

gene_filtering = c("ELOVL6", "TESPA1", "ESR1", "ATF7IP2", "SCAI", "NCOA3", "SLX4IP", "HDAC8", "LCLAT1", "XRCC4",
    "WDR27", "TDRD3", "ZNF678", "PIGN", "ARFIP1", "EYA3", "OPRM1", "FAR2", "SPOPL", "PEX14", "NLK",
    "BTRC", "ZNF407", "NSMCE2", "GLCCI1", "SLC10A7", "FOXN2", "METTL8", "CCDC146", "GRK5", "MBNL1",
    "PPP3CA", "TMEM164","HERPUD2", "MSC−AS1", "IKZF2", "PRKAG2", "PRAG1", "CLOCK", "LRRC8D", "VPS54", "SNX10", "ZDHHC21",
    "BTBD11", "CHD7", "SERINC5", "INPP5A", "RGPD2", "AGK", "WDSUB1", "SORBS1", "CYTH3", "DYNC2H1",
    "ESYT2", "TSPAN18", "NRIP1", "EPB41L2", "PDE9A", "P2RY14", "BACH2", "VWDE", "EMP1", "IL1RAP",
    "BCL2", "FNDC3B", "GNPTAB", "UBE2F", "LIPC", "CHSY3", "ENOX1", "MAP4K3−DT", "ATP9A", "CD70","KLRD1", "KLRC3"
     )
filtered_genes = setdiff(genes, gene_filtering)

df_wide_filtered <- df_wide[filtered_tfs_v2, filtered_genes]
corr_matrix_filtered <- cor(t(df_wide_filtered),method="pearson")

# Create the color function (same as pheatmap default)
pheatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
piyg_colors <- colorRampPalette(rev(brewer.pal(n = 9, name ="PiYG")))(100)

# Define color function with 0 as the center
col_fun <- colorRamp2(
  breaks = c(-1, 0, 1),
  colors = piyg_colors[c(1, ceiling(length(piyg_colors)/2), length(piyg_colors))]
)

#---------------------------------------------------------------------------------------
# Heatmap1
row_clust <- cutree(hclust(dist(corr_matrix_filtered)), k = 4)
col_clust <- cutree(hclust(dist(t(corr_matrix_filtered))), k = 4)

idx_cluster1 <- names(row_clust)[row_clust == 1]

# Subcluster cluster 1 elements
subclust <- cutree(hclust(dist(corr_matrix_filtered[idx_cluster1, ])), k = 2)

# Rebuild cluster assignment
row_split <- row_clust
row_split[idx_cluster1] <- paste0("1", letters[subclust])  # "1a", "1b"
row_split <- factor(row_split)  # ensure it's a factor

ht <- ComplexHeatmap::Heatmap(matrix = corr_matrix_filtered,
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              use_raster = TRUE,
              row_names_gp = grid::gpar(fontsize = 2),
              column_names_gp = grid::gpar(fontsize = 2),
              row_split = row_split,
              column_split = row_split
              )

pdf('Reproducibility/Results/Plots/CD8_T_NK_ILC/Figure5K_ht1.pdf', width = 8, height = 8)
 ht
dev.off()

#---------------------------------------------------------------------------------------
# Heatmap2

ht_drawn <- draw(ht)
TF_order_list <- row_order(ht_drawn)
flat_vector <- unname(do.call(c, TF_order_list))
ordered_TF_names <- rownames(corr_matrix_filtered)[flat_vector]

corr_matrix_filtered_inv <- cor(df_wide_filtered,method="pearson")

# Get original clustering with k = 3
row_clust <- cutree(hclust(dist(corr_matrix_filtered_inv)), k = 3)
col_clust <- cutree(hclust(dist(t(corr_matrix_filtered_inv))), k = 3)

idx_cluster1 <- names(row_clust)[row_clust == 1]

# Subcluster cluster 1 elements
subclust <- cutree(hclust(dist(corr_matrix_filtered_inv[idx_cluster1, ])), k = 2)

# Rebuild cluster assignment
row_split2 <- row_clust
row_split2[idx_cluster1] <- paste0("1", letters[subclust])
row_split2 <- factor(row_split2)

ht2 = ComplexHeatmap::Heatmap(corr_matrix_filtered_inv,
         col = pheatmap_colors,
         cluster_rows = FALSE,
         cluster_columns = FALSE,
         use_raster = TRUE,
         row_names_gp = grid::gpar(fontsize = 1),
         column_names_gp = grid::gpar(fontsize = 1),
         row_split = row_split2,
         column_split = row_split2
         )

pdf('Reproducibility/Results/Plots/CD8_T_NK_ILC/Figure5K_ht2.pdf', width = 8, height = 8)
 ht2
dev.off()

#---------------------------------------------------------------------------------------
# Heatmap3

ht2_drawn <- draw(ht2)
gene_order_list <- row_order(ht2_drawn)
flat_vector2 <- unname(do.call(c, gene_order_list))
ordered_gene_names <- rownames(corr_matrix_filtered_inv)[flat_vector2]

df_wide_filtered2 = df_wide_filtered
df_wide_filtered2[df_wide_filtered2 <0] =0

column_split <- rep(c("1a", "1b", "2", "3"), times = c(220, 98, 249, 222))
column_split <- factor(column_split, levels = c("1a", "1b", "2", "3"))
names(column_split) <- ordered_gene_names

row_split <- rep(c("1a", '1b', "2", "3", "4"), times = c(18, 27, 17, 23, 33))
row_split <- factor(row_split, levels = c("1a", '1b', "2", "3", "4"))
names(row_split) <- ordered_TF_names

ht3 = ComplexHeatmap::Heatmap(df_wide_filtered2[ordered_TF_names, ordered_gene_names],
         col = viridis(100),
         cluster_rows = FALSE,
         cluster_columns = FALSE,
         use_raster = TRUE,
         row_names_gp = grid::gpar(fontsize = 2),
         column_names_gp = grid::gpar(fontsize = 2),
         column_split = column_split,
         row_split = row_split
         )

pdf('Reproducibility/Results/Plots/CD8_T_NK_ILC/Figure5K_ht3.pdf', width = 8, height = 8)
 ht3
dev.off()


#########################################################
## Fig.5L  TF-gene program alluvium
#########################################################

Gene_1a = ordered_gene_names[1:220]
Gene_1b = ordered_gene_names[221:318]
Gene_2 = ordered_gene_names[319:567]
Gene_3 = ordered_gene_names[568:789]

TF_1a = ordered_TF_names[1:18]
TF_1b = ordered_TF_names[19:45]
TF_2 = ordered_TF_names[46:62]
TF_3 = ordered_TF_names[63:85]
TF_4 = ordered_TF_names[86:118]

gene_modules <- setNames(
  c(
    rep("Gene_1a", 220),
    rep("Gene_1b", 98),
    rep("Gene_2", 249),
    rep("Gene_3", 222)
  ),
  ordered_gene_names[1:789]
)

tf_modules <- setNames(
  c(
    rep("TF_1a", 18),
    rep("TF_1b", 27),
    rep("TF_2", 17),
    rep("TF_3", 23),
    rep("TF_4", 33)
  ),
  ordered_TF_names[1:118]
)

importance_matrix = df_wide_filtered %>% t()

importance_long <- as.data.frame(as.table(importance_matrix))
colnames(importance_long) <- c("Gene", "TF", "Importance")

importance_long$Gene_Module <- gene_modules[as.character(importance_long$Gene)]
importance_long$TF_Module <- tf_modules[as.character(importance_long$TF)]

# Aggregate regulation strength
module_regulation <- importance_long %>%
  dplyr::filter(!is.na(Gene_Module), !is.na(TF_Module)) %>%
  group_by(TF_Module, Gene_Module) %>%
  summarise(
    Strength = sum(Importance) / (n_distinct(TF) * n_distinct(Gene)),
    .groups = "drop"
  )

gg_df <- module_regulation %>%
  dplyr::filter(Strength > 0.1) %>%
  dplyr::rename(axis1 = TF_Module,
         axis2 = Gene_Module,
         weight = Strength)

p = ggplot(gg_df,
       aes(axis1 = axis1, axis2 = axis2, y = weight)) +
  geom_alluvium(aes(fill = axis1), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("TF Module", "Gene Module"), expand = c(.05, .05)) +
  scale_fill_manual(values = c('#D81B60', '#1E88E5', '#FFC107', '#004D40', '#F4511E')) +
  labs(title = "TF Module → Gene Module Regulation Strength",
       y = "Strength") +
  theme_minimal()

pdf('Reproducibility/Results/Plots/CD8_T_NK_ILC/Figure5L.pdf', width = 8, height = 8)
 p
dev.off()


#########################################################
## Fig.S11G Violin plot of CD8 GPs in NK cells
#########################################################
sig_df = fread_n('Reproducibility/Results/VISIONR/UC_DOGMA_CD8_T_NK_ILC_signature_score_scenicplue_module.txt')
sig_df = sig_df[colnames(DOGMA),]
NK_celltypes = c("NK_CD56_CD49a_Hi_CD103_Hi","NK_CD56_CD49a_Hi_CD103_Lo","NK_CD56_CD49a_Lo","NK_CD56_dim")

Idents(DOGMA) = 'celltype'
DOGMA[["signature"]] = CreateAssayObject(data = sig_df %>% t() %>% as.matrix())
DOGMA_tmp = subset(DOGMA, ident = NK_celltypes)

DefaultAssay(DOGMA_tmp) <- "signature"
sig_mat <- DOGMA_tmp[["signature"]]@data
genes   <- rownames(sig_mat)
groups  <- DOGMA_tmp@meta.data$celltype

# Build a long dataframe: one row per (cell, gene)
df_long <- as.data.frame(t(as.matrix(sig_mat))) %>%
  tibble::rownames_to_column(var = "cell") %>%
  mutate(Group = groups) %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "Gene",
    values_to = "Value"
  ) %>%
  # keep only these two groups, and order them
  dplyr::filter(Group %in% NK_celltypes) %>%
  mutate(
    Group = factor(Group, levels = NK_celltypes),
    Gene  = factor(Gene, levels = genes)  # facet order = rownames(sig assay)
  )

# Colors
colors <- c("#42a580","#ac484d","#6a8db4","#cb4b7c")

# Plot (one figure with facets)
p <- ggplot(df_long, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  labs(x = NULL, y = "Signature score") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~ Gene, scales = "free_y") 

pdf_3('Reproducibility/Results/Plots/CD8_T_NK_ILC/FigureS11G.pdf', h = 4, w = 7)
 plot(p)
dev.off()


#########################################################
## Fig.S11H Box plot of NK cells fraction
#########################################################

BCG_pt = c('BC_048','BC_039','BC_027','BC_044','BC_043','BC_032','BC_040','BC_047')
Normal_pt = c('c_cell-p2_N-bulk','c_cell-p3_N-bulk','c_cell-p4_N-bulk','c_cell-p5_N-bulk','nat_com-p1_N-bulk','nat_com-p3_N-bulk','nat_com-p6_N-bulk')

df_NK = fread_n("Reproducibility/Results/scANVI/BC/TNK/Atlas_level_integration_TNK_metadata.txt") %>%
        mutate(celltype = if_else(paper != "DOGMA", predicted_celltype, celltype)) %>%
        mutate(patient = if_else(paper != "DOGMA", batch_id, sample)) %>%
        dplyr::filter(., celltype %in% c('NK_CD56_CD49a_Hi_CD103_Hi','NK_CD56_CD49a_Hi_CD103_Lo','NK_CD56_CD49a_Lo','NK_CD56_dim')) %>%
        dplyr::filter(., !patient %in% BCG_pt)

# Filter the patient n < 50
df_NK_filtered <- df_NK %>%
  group_by(patient) %>%
  dplyr::filter(n() >= 50) %>%
  ungroup() %>% as.data.frame()

df_NK$celltype2 = df_NK$celltype %>% 
                  fct_collapse(., NK_CD56_CD49a_Hi = c('NK_CD56_CD49a_Hi_CD103_Hi','NK_CD56_CD49a_Hi_CD103_Lo'))

df_NK$CB = rownames(df_NK)
df_NK$patient = df_NK$patient[,drop=TRUE]
cell_number_df = table(df_NK[,"patient"], df_NK[,"celltype2"]) %>% as.data.frame.matrix(.)

prop_df = data.frame(100*cell_number_df/rowSums(cell_number_df)) %>% 
                     rownames_to_column("patient") %>% 
                     pivot_longer(-patient, names_to = "celltype", values_to = "prop")

prop_df <- prop_df %>%
  dplyr::mutate(primary = ifelse(patient %in% Normal_pt, 'Normal', 'Primary'))

p <- ggplot(prop_df, aes(x = primary, y = prop)) + 
  facet_wrap(~celltype) +
  stat_boxplot(geom = "errorbar", width = 0.4) + 
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  geom_jitter(aes(fill = primary), shape = 21, color = "black", size = 4) + # Color outline black, fill based on primary
  scale_fill_manual(values = c("#0073C299", "#EFC00099")) + # Set custom fill colors
  labs(x = "Stage", y = "Value") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 20)) + 
  ylab("% Total NK cells") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(breaks = seq(0, max(prop_df$prop, na.rm = TRUE), by = 25)) + # Set y-axis breaks by 25
  stat_compare_means(aes(group = primary),
                     comparisons = list(c("Normal", "Primary")),
                     method = "wilcox.test",
                     label = "p.signif",
                     label.y = 105, # Adjust as needed
                     hide.ns = FALSE,
                     tip.length = 0)

paste0('Reproducibility/Results/Plots/CD8_T_NK_ILC/FigureS11H.pdf') %>% pdf_3(., w=12, h=5)
 plot(p)
dev.off()
