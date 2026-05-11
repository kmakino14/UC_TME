#####################################################
# DOGMA/Slide-seqV2 analysis - ECOTYPE -
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
suppressMessages(library(multinichenetr))
suppressMessages(library(CellChat))

suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(ggrepel))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(ggpubr))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(SCpubr))
suppressMessages(library(corrplot))

suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsci))
suppressMessages(library(colorspace))
suppressMessages(library(scales))

suppressMessages(library(escape))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(escape))
suppressMessages(library(forestmodel))

set.seed(1234)

Blue_Red = rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7", "#FFFFFF", "#D1E5F0",
                 "#92C5DE", "#4393C3", "#2166AC", "#053061"))
Orange_Purple = rev(brewer.pal(11, "PuOr"))

#################################################
# Fig.7A  ECOTYPE DOGMA-seq data
#################################################
metadata = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_metadata.txt")
FACS_prop = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_FACS_and_clinical_info.txt")

df_immune = dplyr::filter(metadata, lineage %in% c('B','CD4_T','CD8_T_NK_ILC','Myeloid'))
df_non_immune = dplyr::filter(metadata, !lineage %in% c('B','CD4_T','CD8_T_NK_ILC','Myeloid'))

res_immune = as.data.frame.matrix(table(df_immune$sample, df_immune$celltype)) * FACS_prop$CD45_pos
res_non_immune = as.data.frame.matrix(table(df_non_immune$sample, df_non_immune$celltype)) * FACS_prop$CD45_neg

res <- cbind(res_immune, res_non_immune) %>% as.data.frame()
df_prop <- res %>%
  { 100 * . / rowSums(.) } %>%
  as.data.frame()

df_prop['Epithelial'] = df_prop['LUM'] + df_prop["MES"] + df_prop["NEC"] + df_prop["Normal"] + df_prop["NRP"] + df_prop["SQM"]
df_prop = df_prop %>%
  mutate(
    iCAF   = iCAF_SLC14A1 + iCAF_CD321 + proCAF,
    matCAF = matCAF + matCAP,
    CD8_Tex = CD8_Tex_1 + CD8_Tex_2,
    NK_CD56_CD49a_Hi = NK_CD56_CD49a_Hi_CD103_Hi + NK_CD56_CD49a_Hi_CD103_Lo

  ) %>%
  dplyr::select(-any_of(c("iCAF_SLC14A1","iCAF_CD321","matCAP","CD8_Tex_1","CD8_Tex_2",'LUM','MES','NEC','Normal','NRP','SQM',
                          "NK_CD56_CD49a_Hi_CD103_Hi","NK_CD56_CD49a_Hi_CD103_Lo",'proCAF','Mast'))) %>%
  t()

df_norm = rank_norm(df_prop)
df_norm2 <- df_norm[!rownames(df_norm) %in% "Epithelial", ]

pheatmap_result = plot_corrmatrix(df = df_norm2, col_list = Blue_Red,
                cor_method = 'spearman',
                save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure7A_celltype_module.pdf")

#################################################################
df_cor = cor(df_norm2 %>% as.matrix() %>% t() ,method='spearman')

row.order = c('Epithelial', rownames(df_norm2[pheatmap_result$tree_row$order,]))
clustered_data <- df_norm[row.order,]

col_anno_df = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_FACS_and_clinical_info.txt") %>% as.data.frame() %>%
              dplyr::select(., c("status","main_histology"))
col_anno_df$status = factor(col_anno_df$status, levels = c("Early_L","Early_H","Early_UTUC","Advanced_BC","Advanced_UTUC","post_BCG"))
col_anno_df$main_histology = factor(col_anno_df$main_histology,levels = c("UC","Small_cell","Inverted_papilloma","No_malignancy"))
row_anno_df = data.frame(ECOTYPE = rep(c('Tumor','EC1','EC2','EC3','EC4','EC5','EC6','EC7','EC8','EC9', 'EC10'), 
                         times = c(1,4,8,3,7,4,2,7,2,3,2)))
rownames(row_anno_df) = rownames(clustered_data)

# color 
pal1 <- ggsci::scale_color_d3('category10')
pal2 <- ggsci::scale_color_tron()

col.group <- pal1$palette(6)
col.histo <- pal2$palette(4)
col.eco <- hue_pal()(12)

names(col.group) <- c("Early_L","Early_H","Early_UTUC","Advanced_BC","Advanced_UTUC","post_BCG")
names(col.histo) <- c("UC","Small_cell","Inverted_papilloma","No_malignancy")
names(col.eco) <- c('Tumor','EC1','EC2','EC3','EC4','EC5','EC6',"EC7","EC8",'EC9','EC10','EC11')

anno_color <- list(
  status = col.group, 
  main_histology = col.histo,
  ECOTYPE = col.eco
  )

# Heatmap
heatmap_dendro = pheatmap::pheatmap(clustered_data,
                                    annotation_row = row_anno_df,
                                    annotation_col = col_anno_df,
                                    annotation_colors = anno_color, 
                                    cluster_rows = FALSE,
                                    cluster_cols = TRUE,
                                    cellwidth = 10, cellheight = 10, color=cividis(50), 
                                    clustering_method = "ward.D2", breaks = seq(-1,2,length.out=50),
                                    filename=paste0("Reproducibility/Results/Plots/ECOTYPE/Figure7A_celltype_abundance_concat.pdf"), 
                                    border_color = "black" ,na_col = "gray90",
                                    main = paste0("ranknorm_spearman_clustered_ecotype_w_dendrogram"))

after_clustering_order_Pts <- heatmap_dendro$tree_col$order

pheatmap::pheatmap(clustered_data[,after_clustering_order_Pts],
                   annotation_row = row_anno_df,
                   annotation_col = col_anno_df,
                   annotation_colors = anno_color, 
                   gaps_row=c(1,5,13,16,19,23,27,29,36,38,41),
                   gaps_col=c(8,16,21,30,34),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   cellwidth = 10, cellheight = 10, color=cividis(50), 
                   clustering_method = "ward.D2", breaks = seq(-1,2,length.out=50),
                   filename=paste0("Reproducibility/Results/Plots/ECOTYPE/Figure7A_celltype_abundance_concat_annotate.pdf"), 
                   border_color = "black" ,na_col = "gray90",
                   main = paste0("ranknorm_spearman_clustered_ecotype"))


##################################################
# Fig.S13E  ECOTYPE TNK external validation
##################################################
########################
## External dataset
########################
lineage_use <- c("TNK")

integrated_meta = fread_n("Reproducibility/Results/scANVI/BC/TNK/Atlas_level_integration_TNK_metadata.txt") %>% 
                  dplyr::filter(., !paper %in% c('DOGMA')) %>%
                  dplyr::filter(., predicted_lineage %in% lineage_use)

res = table(integrated_meta[,"batch_id"],integrated_meta[,"predicted_celltype"]) %>% as.data.frame.matrix()
res$CD8_T_proliferative = NULL
res$CD4_T_proliferative = NULL
res$CD8_Tm = res$CD8_Tcm + res$CD8_Trm
res$CD8_Tcm = NULL
res$CD8_Trm = NULL

res$sum = rowSums(res)
res = dplyr::filter(res, sum > 1500) # 2000 
res$sum = NULL

df_prop = data.frame(100*res/rowSums(res)) %>% as.data.frame() %>% t()  # celltype*pts
df_ranknorm = rank_norm(df_prop)

plot_corrmatrix(df = df_ranknorm, col_list =Orange_Purple,
                cor_method = 'spearman',
                clustering_method = 'ward.D2',
                save_path ="Reproducibility/Results/Plots/ECOTYPE/FigureS13E_external.pdf")

########################
## In house dataset
########################
metadata = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_metadata.txt")
df_TNK = dplyr::filter(metadata, coarse_celltype %in% c('CD4_Tconv','Treg',"CD8_T",'NK_ILC')) %>%
         dplyr::filter(!celltype %in% c('CD4_T_proliferative','CD8_T_proliferative'))

res = table(df_TNK[,"sample"],df_TNK[,"celltype"]) %>% as.data.frame.matrix()
res$CD8_Tm = res$CD8_Tcm + res$CD8_Trm
res$CD8_Tcm = NULL
res$CD8_Trm = NULL
df_prop = data.frame(100*res/rowSums(res)) %>% as.data.frame() %>% t() # celltype*pts
df_ranknorm = rank_norm(df_prop)

plot_corrmatrix(df = df_ranknorm, col_list = Blue_Red,
                cor_method = 'spearman',
                clustering_method = 'ward.D2',
                save_path ="Reproducibility/Results/Plots/ECOTYPE/FigureS13E_in_house.pdf")


#################################################
# Fig.7B/Fig.S13F  Signature score boxplot
#################################################
# Load 
sig_df_1 = fread_n("Reproducibility/Results/VISIONR/UC_DOGMA_Malignant_signature_score_literature.txt")
sig_df_2 = fread_n("Reproducibility/Results/VISIONR/UC_DOGMA_Malignant_signature_score_hotspot.txt")
sig_df = cbind(sig_df_1,sig_df_2)
metadata = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_metadata.txt") %>% dplyr::filter(., celltype %in% c("LUM",'NRP','SQM','MES','NEC'))
metadata = metadata[rownames(sig_df),]

sig_df_w_meta = cbind(sig_df, metadata) %>% dplyr::select(.,c(colnames(sig_df), 'sample', 'group')) 
sig_df_avg <- sig_df_w_meta %>%
  group_by(sample) %>%
  summarise(across(-group, mean, na.rm = TRUE)) %>% as.data.frame()

sig_df_avg=sig_df_avg %>% column_to_rownames('sample')
sig_df_norm <- as.data.frame(apply(sig_df_avg, 2, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

plot_tme_boxplots <- function(df, alt, tmp_feature, tmp_TME, save_path){
  tme_sets <- list(TME_1 = TME_1, TME_3 = TME_3, TME_4 = TME_4, TME_5 = TME_5, TME_6 = TME_6)
  colors <- c("TME_1"="#bd9e39","TME_3"="#9ecae1","TME_4"="#d6616b",
              "TME_5"="#a1d99b","TME_6"="#ce6dbd","Others"="#DAD8C9")
  current <- tme_sets[[tmp_TME]]
  others  <- setdiff(unique(unlist(tme_sets)), current)

  # Assign group labels
  df$Group <- NA_character_
  df$Group[rownames(df) %in% current] <- tmp_TME
  df$Group[rownames(df) %in% others]  <- "Others"

  # Long format
  df_long <- df %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(cols = -c(Sample, Group),
                        names_to = "Feature", values_to = "Value") %>%
    dplyr::filter(!is.na(Group)) %>%
    dplyr::filter(Feature %in% tmp_feature)

  df_long$Group = factor(df_long$Group, levels = c(tmp_TME,'Others'))

  # Plot
  p <- ggplot(df_long, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.6, shape = 16) +
    facet_wrap(~ Feature, scales = "free_y", ncol = 6) +
    stat_compare_means(
      comparisons = list(c(tmp_TME, "Others")),
      method = "wilcox.test",
      method.args = if (is.null(alt)) list() else list(alternative = alt),
      label = "p.signif",
      size = 2.3,
      label.y.npc = 0.90,
      hide.ns = TRUE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    scale_fill_manual(values = colors[c(tmp_TME, "Others")]) +
    theme_classic(base_size = 8) +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "", y = "Score", title = paste0("Wilcoxon: ", tmp_TME, " vs Others"))

  # Save
  pdf(save_path, width = 2, height = 2)
  plot(p)
  dev.off()
}

plot_tme_boxplots_3group <- function(df, alt, tmp_feature, save_path){
  tme_sets <- list(TME_1 = TME_1, TME_3 = TME_3, TME_4 = TME_4, TME_5 = TME_5, TME_6 = TME_6)
  colors <- c("TME_1"="#bd9e39","TME_3"="#9ecae1","TME_4"="#d6616b",
              "TME_5"="#a1d99b","TME_6"="#ce6dbd","Others"="#DAD8C9")

  df$Group <- NA_character_
  df$Group[rownames(df) %in% TME_4] <- "TME_4"
  df$Group[rownames(df) %in% TME_5] <- "TME_5"
  df$Group[!rownames(df) %in% c(TME_4,TME_5)]  <- "Others"

  # Long format
  df_long <- df %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(cols = -c(Sample, Group),
                        names_to = "Feature", values_to = "Value") %>%
    dplyr::filter(!is.na(Group)) %>%
    dplyr::filter(Feature %in% tmp_feature)

  df_long$Group = factor(df_long$Group, levels = c("TME_4","TME_5",'Others'))

  # Plot
  p <- ggplot(df_long, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.6, shape = 16) +
    facet_wrap(~ Feature, scales = "free_y", ncol = 6) +
    stat_compare_means(
      comparisons = list(c('TME_4', "Others"),c('TME_5', "Others"),c('TME_4', "TME_5")),
      method = "wilcox.test",
      method.args = if (is.null(alt)) list() else list(alternative = alt),
      label = "p.signif",
      size = 2.3,
      label.y.npc = 0.90,
      p.adjust.method = "BH",
      hide.ns = TRUE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    scale_fill_manual(values = colors[c('TME_4','TME_5', "Others")]) +
    theme_classic(base_size = 8) +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(x = "", y = "Score", title = paste0("Wilcoxon: TME_4 vs TME_5 vs Others"))

  # Save
  pdf(save_path, width = 3, height = 2)
  plot(p)
  dev.off()
}

plot_tme_boxplots(df = sig_df_norm, alt = NULL, tmp_feature = "Oxidative_phosphorylation", 
                  tmp_TME = "TME_1", save_path = "Reproducibility/Results/Plots/ECOTYPE/FigureS13F_TME_1_vs_others.pdf")
plot_tme_boxplots(df = sig_df_norm, alt = NULL, tmp_feature = "Cytokines_and_inflammatory_response", 
                  tmp_TME = "TME_3", save_path = "Reproducibility/Results/Plots/ECOTYPE/FigureS13F_TME_3_vs_others.pdf")

plot_tme_boxplots_3group(df = sig_df_norm, alt = 'less', tmp_feature = "Luminal", 
                         save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure7B_LUM.pdf")
plot_tme_boxplots_3group(df = sig_df_norm, alt = 'greater', tmp_feature = "Partial_epithelial-mesenchymal_transition", 
                         save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure7B_pEMT.pdf")


#################################################
# Fig.7B  TF activity boxplot
#################################################
TF_activity = fread_n("LINGER/Primary/Malignant/output/cell_population_TF_activity_zscore.txt")

metadata = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_metadata.txt") %>% dplyr::filter(., celltype %in% c("LUM",'NRP','SQM','MES','NEC'))
metadata = metadata[rownames(TF_activity),]

TF_df_w_meta = cbind(TF_activity, metadata) %>% 
               dplyr::select(.,c(colnames(TF_activity), 'sample', 'group')) 
TF_df_avg <- TF_df_w_meta %>%
  group_by(sample) %>%
  summarise(across(-group, mean, na.rm = TRUE)) %>% as.data.frame()

TF_df_avg=TF_df_avg %>% column_to_rownames('sample')
TF_df_norm <- as.data.frame(apply(TF_df_avg, 2, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

TF_df_norm = TF_df_norm[c(TME_1,TME_3,TME_4,TME_5,TME_6), ]

plot_tme_boxplots(df = TF_df_norm, alt = NULL, tmp_feature = "NR2F2", 
                  tmp_TME = "TME_6", save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure7B_NR2F2.pdf")
plot_tme_boxplots_3group(df = TF_df_norm, alt = 'less', tmp_feature = "GATA3", 
                         save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure7B_GATA3.pdf")

#########################################################
## Fig.6J  Spatial celltype enrichment
#########################################################
# enrichment plot
enr = fread_n("Reproducibility/Results/Slide-seqV2/spatial_coarse_celltype_enrichment_by_regions.txt")

cts_order = c("CD4_Tn", "CD4_Tcm", "CD4_Tsen", "CD4_Tfh-like", "CD4_T_proliferative", "Treg",
  "CD8_Tn", "CD8_Tcm", "CD8_Tem", "NK_ILC", "CD4_Th17", "CD8_Temra",
  "B", "CD4_T_CD26",
  "CD8_T_proliferative", "CD4_Tph-like", "Plasma", "CD8_Tex", "CD4_CTL", "CD8_Trm",
  "MoMac", "DC", "MSC")

p = plot_significances(
    significances = enr,
    p_key='p_mwu_fdr_bh', 
    value_key='predictions_filt', 
    group_key='regions', 
    enrichment_key = "enrichment",
    enriched_label = "enriched",
    pmax = 0.01,
    pmin = 1e-5,
    annotate_pvalues = FALSE,
    value_cluster = FALSE,
    value_order = cts_order,
    group_cluster = FALSE,
    group_order = c('Region 1','Region 2','Region 3','Region 4'),
    method = 'ward.D2'
    )

paste0("Reproducibility/Results/Plots/Slide-seqV2/Figure6J_enrichment.pdf") %>% pdf(w=8, h=10)
 p
dev.off()

#################################################
# Fig.S13D  ECOTYPE Slide-seqV2
#################################################
# 1) celltype
cell_counts = fread_n("Reproducibility/Results/Slide-seqV2/spatial_celltype_compostion_by_patient.txt")
df_prop = data.frame(100*cell_counts/rowSums(cell_counts)) %>% as.data.frame() %>% t()  # celltype*pts
df_prop = t(df_prop) %>% as.data.frame()
df_prop$iCAF = df_prop$iCAF + df_prop$proCAF
df_prop$NK_CD56_CD49a_Hi = df_prop$NK_CD56_CD49a_Hi_CD103_Hi + df_prop$NK_CD56_CD49a_Hi_CD103_Lo
df_prop$proCAF = NULL
df_prop$NK_CD56_CD49a_Hi_CD103_Hi = NULL
df_prop$NK_CD56_CD49a_Hi_CD103_Lo = NULL
df_prop = t(df_prop) %>% as.data.frame()

# Replacement dictionary
repl <- c(
  "LUM"          = "Epithelial",
  "CD4_Tfh.like" = "CD4_Tfh-like",
  "CD4_Tph.like" = "CD4_Tph-like",
  "MDSC.like"    = "MDSC-like"
)

# Apply replacements
rownames(df_prop) <- plyr::mapvalues(
  rownames(df_prop),
  from = names(repl),
  to   = repl,
  warn_missing = FALSE
)

df_ranknorm = rank_norm(df_prop)

#################################################
# celltype abundance
EC0 = c('Epithelial')
EC1 = c("GC_B","Atypical_B","B_memory","CD4_Tfh-like")
EC2 = c("pDC","ILC3","NK_CD56_CD49a_Lo","CD4_Tcm","CD4_Tn","Treg_naive","B_naive","CD8_Tn")
EC3 = c("CD8_Tcm","CD4_Tsen","TAM_FOLR2")
EC4 = c("Mono","CD8_Temra","NK_CD56_dim")
EC5 = c("cDC1","mregDC","cDC2","preDC")
EC6 = c("vSMC","contCAP","Endothelial","iCAF")
EC7 = c("TAM_TREM2","matCAF")
EC8 = c("CD4_Tph-like","CD8_Tex","Treg_effector","NK_CD56_CD49a_Hi", "CD8_Trm","CD8_T_proliferative","MDSC-like")
EC9 = c("CD8_Tem","Plasma")
EC10 = c("MAIT","CD4_T_proliferative","CD4_Th17")
EC11 = c("CD4_CTL","CD4_T_CD26")

df_ranknorm = t(df_ranknorm) %>% as.data.frame()
df_ranknorm$Mast <- NULL
df_ranknorm = t(df_ranknorm) %>% as.data.frame()

row.order = c(EC0,EC1,EC2,EC3,EC4,EC5,EC6,EC7,EC8,EC9,EC10,EC11)
clustered_data <- df_ranknorm[row.order,]

EC_list <- list(
  EC0 = EC0, EC1 = EC1, EC2 = EC2, EC3 = EC3, EC4 = EC4, EC5 = EC5,
  EC6 = EC6, EC7 = EC7, EC8 = EC8, EC9 = EC9, EC10 = EC10, EC11 = EC11
)

compute_EC_means <- function(clustered_data, EC_list) {
  # Ensure matrix
  clustered_data <- as.matrix(clustered_data)

  # Compute mean for each EC group
  EC_means <- lapply(names(EC_list), function(ec) {
    rows <- intersect(EC_list[[ec]], rownames(clustered_data))
    if (length(rows) == 0) return(rep(NA, ncol(clustered_data)))  # avoid NULL
    colMeans(clustered_data[rows, , drop = FALSE])
  })

  # Convert list → matrix
  EC_mean_mat <- do.call(rbind, EC_means)
  rownames(EC_mean_mat) <- names(EC_list)

  EC_mean_mat = as.matrix(EC_mean_mat)
  return(EC_mean_mat)
}

EC_mean_mat <- compute_EC_means(clustered_data, EC_list)

# Heatmap
pheatmap::pheatmap(EC_mean_mat,
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   cellwidth = 10, cellheight = 10, color=cividis(50), 
                   clustering_method = "ward.D2", breaks = seq(-1,1.5,length.out=50),
                   filename=paste0("Reproducibility/Results/Plots/Slide-seqV2/FigureS13D.pdf"), 
                   border_color = "black" ,na_col = "gray90",
                   main = paste0("ranknorm_spearman_slide-seqV2"))

############################################################################
# 2) Assign TME subtype
df_sc = fread_n("ECOTYPE/ranknorm_cell_abundance_w_epithelial.txt")
df_sc = t(df_sc) %>% as.data.frame()
df_sc$Mast <- NULL
df_sc = t(df_sc) %>% as.data.frame()

EC_mat_sc <- compute_EC_means(df_sc, EC_list)
EC_mat_sp <- compute_EC_means(clustered_data, EC_list)

## -----------------------------
## 1. Hierarchical clustering on EC_mat_sc (samples)
##    (Skip this part if you already have hc_sc and clusters_sc)
## -----------------------------
# Cut into 6 clusters
clusters_sc <- c(
  # Cluster 6
  BC_013 = 6, BC_011 = 6, BC_045 = 6, BC_038 = 6,
  BC_014 = 6, BC_010 = 6, UTUC_001 = 6, UTUC_007 = 6,

  # Cluster 5
  BC_050 = 5, BC_026 = 5, BC_022 = 5, BC_008 = 5,

  # Cluster 4
  BC_042 = 4, BC_018 = 4, BC_041 = 4, BC_028 = 4,
  BC_031 = 4, BC_020 = 4, BC_030 = 4, UTUC_003 = 4, BC_036 = 4,

  # Cluster 1
  BC_029 = 1, UTUC_006 = 1, BC_017 = 1, BC_021 = 1,
  BC_012 = 1, UTUC_004 = 1, BC_019 = 1, BC_024 = 1,

  # Cluster 3
  BC_037 = 3, BC_016 = 3, BC_023 = 3, BC_033 = 3, UTUC_005 = 3,

  # Cluster 2
  BC_047 = 2, BC_039 = 2, BC_044 = 2, BC_032 = 2,
  BC_048 = 2, BC_043 = 2, BC_027 = 2, BC_040 = 2
)

table(clusters_sc)

## -----------------------------
## 2. Compute cluster centroids in ECOTYPE space
##    (mean ECOTYPE profile for each of the 6 TME clusters)
## -----------------------------
cluster_ids <- sort(unique(clusters_sc))

centers_sc <- sapply(cluster_ids, function(cl){
  samp_in_cl <- names(clusters_sc)[clusters_sc == cl]
  rowMeans(EC_mat_sc[, samp_in_cl, drop = FALSE])
})

# centers_sc: matrix with rows = ECOTYPE, cols = clusters
colnames(centers_sc) <- paste0("Cluster", cluster_ids)

## -----------------------------
## 3. Align ECOTYPEs between sc and spatial matrices
## -----------------------------
common_ec <- intersect(rownames(EC_mat_sc), rownames(EC_mat_sp))

EC_sc_sub   <- EC_mat_sc[common_ec, , drop = FALSE]
EC_sp_sub   <- EC_mat_sp[common_ec, , drop = FALSE]
centers_sub <- centers_sc[common_ec, , drop = FALSE]

# Z-score per ECOTYPE using EC_mat_sc
mu  <- rowMeans(EC_sc_sub)
sdv <- apply(EC_sc_sub, 1, sd)
sdv[sdv == 0] <- 1  # avoid division by zero

scale_fun <- function(mat){
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}

EC_sc_scaled   <- scale_fun(EC_sc_sub)
EC_sp_scaled   <- scale_fun(EC_sp_sub)
centers_scaled <- scale_fun(centers_sub)

## -----------------------------
## 4. Assign each spatial sample to the nearest sc cluster
## -----------------------------
# Compute squared Euclidean distance between each spatial sample
# and each cluster center (in scaled ECOTYPE space)
dist_to_centers <- sapply(1:ncol(EC_sp_scaled), function(j){
  sapply(1:ncol(centers_scaled), function(k){
    sum((EC_sp_scaled[, j] - centers_scaled[, k])^2)
  })
})

rownames(dist_to_centers) <- colnames(centers_scaled)  # Cluster1..Cluster6
colnames(dist_to_centers) <- colnames(EC_sp_scaled)    # spatial samples

# For each spatial sample, pick the closest cluster
cluster_sp_idx <- apply(dist_to_centers, 2, which.min)

# Map back to cluster IDs used in EC_mat_sc
cluster_sp <- cluster_ids[cluster_sp_idx]
names(cluster_sp) <- colnames(EC_sp_scaled)

print(cluster_sp)
BC084 BC093 BC101 BC102 BC094 BC073 BC071 BC066 BC096 BC098 
    1     5     5     1     3     6     4     4     2     2 

#########################################################
## Fig.7D  CCI analysis (CellChat)
#########################################################
scRNA = paste0("Reproducibility/Data/Seurat/UC_DOGMA_seurat_obj_Global_RNA_HQ.rds") %>% readRDS()

scRNA$coarse_celltype2 = fct_collapse(scRNA$coarse_celltype, 
    Tumor = c('UC','NEC'),  
    MSC = c('CAF','Pericyte')
)

scRNA$samples = scRNA$sample
scRNA$samples = factor(scRNA$samples)
scRNA$group = fct_collapse(scRNA$sample,
                TME4 = c('BC_018',"BC_041","BC_028","BC_031","BC_020","BC_030","UTUC_003","BC_036","BC_042"),
                TME6 = c("BC_010","BC_014","BC_038","BC_045","BC_011","BC_013","UTUC_001","UTUC_007"),
                TME5 = c("BC_008","BC_026","BC_022","BC_050"),
                TME3 = c("UTUC_005","BC_033","BC_023","BC_016","BC_037"),
                TME2 = c("BC_040","BC_027","BC_043","BC_048","BC_032","BC_044","BC_039","BC_047"),
                TME1 = c("BC_029","UTUC_006","BC_017","BC_021","BC_012","UTUC_004","BC_019","BC_024")
                )

###############################################################################
Idents(scRNA) = 'group'

for(tmp_group in c('TME4','TME5','TME6')){
  print(paste0(tmp_group, ' start'))

  scRNA_tmp_group = subset(scRNA, idents = tmp_group)

  # 1) Data input & processing and initialization of CellChat object
  cellchat <- createCellChat(object = scRNA_tmp_group, group.by = "coarse_celltype2", assay = "RNA")
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

  saveRDS(cellchat, file = paste0("Reproducibility/Results/CellChat/UC_DOGMA_",tmp_group,".rds"))

  print(paste0(tmp_group, ' finish'))
}

###############################################################################
netVisual_heatmap_fixed <- function(object, signaling = NULL, slot.name = "netP", color.use = NULL, 
                                    color.heatmap = "Reds", title.name = NULL, font.size = 8, font.size.title = 10, 
                                    cluster.rows = FALSE, cluster.cols = FALSE,
                                    remove.isolate = FALSE, 
                                    max.dataset = NULL, 
                                    max.barplot = NULL){
  
  cells.level <- levels(object@idents)
  prob_all <- slot(object, slot.name)$prob
  
  if (signaling %in% dimnames(prob_all)[[3]]) {
    net.diff <- prob_all[,,signaling]
  } else {
    net.diff <- matrix(0, nrow = length(cells.level), ncol = length(cells.level), 
                       dimnames = list(cells.level, cells.level))
  }
  
  mat <- matrix(0, nrow = length(cells.level), ncol = length(cells.level), dimnames = list(cells.level, cells.level))
  common_rows <- intersect(rownames(net.diff), cells.level)
  common_cols <- intersect(colnames(net.diff), cells.level)
  mat[common_rows, common_cols] <- net.diff[common_rows, common_cols]

  if (is.null(max.dataset) || max.dataset <= 0) { 
    max_val_use <- max(mat, na.rm = TRUE) 
  } else { 
    max_val_use <- max.dataset 
  }
  if (max_val_use <= 0) max_val_use <- 0.01
  
  pal_cols <- brewer.pal(n = 9, name = color.heatmap)
  color.heatmap.use = colorRamp2(seq(0, max_val_use, length.out = length(pal_cols)), pal_cols)

  if (is.null(color.use)) {
    if ("colors" %in% slotNames(object)) { color.use <- object@colors } 
    else { color.use <- hcl(h = seq(15, 375, length = length(cells.level) + 1), l = 65, c = 100)[1:length(cells.level)] }
    names(color.use) <- cells.level
    color.use <- color.use[colnames(mat)]
  }

  if (is.null(max.barplot) || max.barplot <= 0) {
    max_bar_use <- max(c(rowSums(mat), colSums(mat)), na.rm = TRUE)
  } else {
    max_bar_use <- max.barplot
  }
  if (max_bar_use <= 0) max_bar_use <- 0.1

  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(mat), border = FALSE, 
                                              gp = gpar(fill = color.use, col=color.use),
                                              ylim = c(0, max_bar_use)), 
                      show_annotation_name = FALSE)
  
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat), border = FALSE, 
                                                  gp = gpar(fill = color.use, col=color.use),
                                                  ylim = c(0, max_bar_use)), 
                          show_annotation_name = FALSE)

  df_label <- data.frame(group = colnames(mat)); rownames(df_label) <- colnames(mat)
  col_anno <- HeatmapAnnotation(df = df_label, col = list(group = color.use), which = "column", 
                                show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = grid::unit(0.2, "cm"))
  row_anno <- HeatmapAnnotation(df = df_label, col = list(group = color.use), which = "row", 
                                show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = grid::unit(0.2, "cm"))

  mat[mat == 0] <- NA
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Prob.",
                bottom_annotation = col_anno, left_annotation = row_anno, 
                top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows, cluster_columns = cluster.cols, 
                row_names_side = "left", row_names_gp = gpar(fontsize = font.size),
                column_names_gp = gpar(fontsize = font.size),
                column_title = title.name, column_title_gp = gpar(fontsize = font.size.title))
  
  return(ht1)
}

Pathway_communication_prob_heatmap <- function(object.list, pathways.show, output_path){
for (i in 1:length(object.list)) {
  object.list[[i]]@idents <- factor(object.list[[i]]@idents, levels = my_levels)
}

pathways.show <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

pdf(output_path, width = 12, height = 7)

for (p in 1:length(pathways.show)) {
  target_pathway <- pathways.show[p]
  current_vals <- c()
  current_sums <- c()
  
  for (i in 1:length(object.list)) {
    if (target_pathway %in% dimnames(object.list[[i]]@netP$prob)[[3]]) {
      m <- object.list[[i]]@netP$prob[,,target_pathway]
      current_vals <- c(current_vals, as.vector(m))
      current_sums <- c(current_sums, rowSums(m), colSums(m))
    }
  }
  
  page_max_val <- if(length(current_vals) > 0) max(current_vals, na.rm = TRUE) else 0.01
  page_max_sum <- if(length(current_sums) > 0) max(current_sums, na.rm = TRUE) else 0.1
  
  if (page_max_val <= 0) page_max_val <- 0.01
  if (page_max_sum <= 0) page_max_sum <- 0.1

  ht <- vector("list", length(object.list))
  for (i in 1:length(object.list)) {
    ht[[i]] <- netVisual_heatmap_fixed(
      object.list[[i]], 
      signaling = target_pathway, 
      color.heatmap = "YlGnBu",
      max.dataset = page_max_val,
      max.barplot = page_max_sum, 
      title.name = paste(target_pathway, names(object.list)[i]),
      cluster.rows = FALSE, 
      cluster.cols = FALSE
    )
  }

  draw(ht[[1]] + ht[[2]], ht_gap = unit(1, "cm"))
}

dev.off()
}

my_levels <- c('Tumor','MSC','Endothelial',"CD4_Tconv","CD8_T","Treg",
               "NK_ILC","B","Mono_Mac","DC","Mast")

###############################################################################
cellchat.TME4 <- readRDS("Reproducibility/Results/CellChat/UC_DOGMA_TME4.rds")
cellchat.TME5 <- readRDS("Reproducibility/Results/CellChat/UC_DOGMA_TME5.rds")
cellchat.TME6 <- readRDS("Reproducibility/Results/CellChat/UC_DOGMA_TME6.rds")

# TME4 vs TME5
object.list <- list(TME4 = cellchat.TME4, TME5 = cellchat.TME5)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
output_path <- "Reproducibility/Results/Plots/Figure7D_TME4_vs_TME5.pdf"
Pathway_communication_prob_heatmap(object.list, pathways.show = c("TGFb","GALECTIN","Netrin"), output_path)

# TME6 vs TME5
object.list <- list(TME6 = cellchat.TME6, TME5 = cellchat.TME5)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
output_path <- "Reproducibility/Results/Plots/Figure7D_TME6_vs_TME5.pdf"
Pathway_communication_prob_heatmap(object.list, pathways.show = c("TGFb","ANGPTL","ncWNT"), output_path)


#################################################
## Fig.7F  Multi Cox regression (BRSpred cohort) 
#################################################

module_df = fread_n("Reproducibility/Data/Signature_gene_list/TableS3_module_genes_Malignant_hotspot.txt")
module_df$Module = fct_collapse(module_df$Annotation,
                                LUM = c("Luminal"),
                                NRP = c("Neural-like progenitor"),
                                SQM = c("Squamous"),
                                MES = c("Mesenchymal"),
                                NEC = c("Neuroendocrine"),
                                MTN = c("Metanephros"),
                                CYG = c("Cycling/G2M"),
                                IFN = c("Interferon signaling"),
                                SEC = c("Secretion"),
                                STR = c("Stress"),
                                pEMT = c("Partial epithelial-mesenchymal transition")
                                )

load(".../CohortA_pre.RData")
load(".../CohortB.RData")
metadata = fread_n(".../Metadata_from_Table_S2.txt") %>%
           dplyr::filter(., Tumor %in% c("Primary"))

common_genes <- intersect(rownames(CohortA_pre), rownames(CohortB))
CohortA_subset <- CohortA_pre[common_genes, ]
CohortB_subset <- CohortB[common_genes, ]
combined_counts <- cbind(CohortA_subset, CohortB_subset)

na_indices <- which(is.na(combined_counts), arr.ind = TRUE)
res <- data.frame(
  Gene = rownames(combined_counts)[na_indices[, 1]],
  Sample = colnames(combined_counts)[na_indices[, 2]]
)
genes_to_remove <- unique(res$Gene)
combined_counts_clean <- combined_counts[!(rownames(combined_counts) %in% genes_to_remove), ]

###################
Get_signature_score_z = function(module_df){
  counts_tmp = combined_counts_clean
  metadata = metadata

  ## ssGSEA
  counts = counts_tmp
  se <- CreateSeuratObject(counts = counts)
  geneSets <- split(rownames(module_df), module_df$Module)
  df_ssGSEA <- enrichIt(obj = se[["RNA"]]@counts, 
                        gene.sets = geneSets, 
                        groups = 3000, 
                        cores = 2, 
                        min.size = 1, 
                        ssGSEA.norm = TRUE)
  df_ssGSEA_z <- apply(df_ssGSEA, 2, function(x) (x - mean(x)) / sd(x))
  df_ssGSEA_z <- as.data.frame(df_ssGSEA_z)
  df_ssGSEA_w_meta = cbind(df_ssGSEA_z[rownames(metadata),],metadata)
  return(df_ssGSEA_w_meta)
}

BRS_cohort = Get_signature_score_z(module_df = module_df)

BRS_cohort$Sex = fct_collapse(BRS_cohort$Sex,
                 Female = c('f','Female'), Male = c('m', 'Male'))
df_for_cox = BRS_cohort %>% 
             dplyr::select(c("CYG","IFN","LUM","MES","NEC","pEMT","SEC","SQM","STR","NRP","Age",
                           "Sex","LVI","Progression","Time_to_prog_or_FUend")) %>%
             transmute(Time_to_prog_or_FUend,
               Progression,
               Age,
               Sex, LVI, 
               LUM,NRP,SQM,MES,NEC,CYG,IFN,SEC,STR,pEMT
               )
plot = forest_model(coxph(Surv(Time_to_prog_or_FUend, Progression) ~ ., df_for_cox))

paste0("Reproducibility/Results/Plots/ECOTYPE/Figure7F.pdf") %>% pdf_3(., w=10,h=10)
  plot
dev.off()


#########################################################
## Fig.S13K  DPP4 expression external dataset
#########################################################

seurat = readRDS(".../GSE269877_dta_cancer.submission.rds")
meta_npj = fread_n("Reproducibility/Results/scANVI/BCG/annotation_transfer_metadata.txt") 
rownames(meta_npj) <- sub("-1$", "", rownames(meta_npj))
seurat@meta.data = meta_npj

celltype = c("CD4_Tn","CD4_Tcm","CD4_Tsen","CD4_T_CD26","CD4_CTL","CD4_Th17","CD4_Tfh-like",
             "CD4_Tph-like","CD4_T_proliferative","Treg_naive","Treg_effector")

Idents(seurat) = 'predicted_celltype'
seurat = subset(seurat, idents = celltype)
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)

Idents(seurat) = 'predicted_celltype'
seurat$predicted_celltype = factor(seurat$predicted_celltype, levels = rev(celltype))
p <- SCpubr::do_ExpressionHeatmap(sample = seurat,
                                  slot = "scale.data",
                                  features = c('SELL',"DPP4",'CCR6','IFNG'),
                                  features.order = c('SELL',"DPP4",'CCR6','IFNG'),
                                  enforce_symmetry = TRUE,
                                  diverging.palette = "RdBu",
                                  cluster = FALSE)

pdf("Reproducibility/Results/Plots/ECOTYPE/FigureS13K.pdf", w=5, h=5)
 plot(p)
dev.off()

#########################################################
## Fig.7G/S13L  Malignant signature score external 
#########################################################
## pEMT

tmp_path = "Reproducibility/Results/scANVI/BCG/External_dataset_malignant_signature_score_scanpy.txt"
sig_df = fread_n(tmp_path)
meta_npj = fread_n("Reproducibility/Results/scANVI/BCG/annotation_transfer_metadata.txt") %>%
           dplyr::filter(., predicted_celltype %in% c("Epithelial"))

sig_df$sample = meta_npj$sample

# Define groups
Naive_w = c("BC_19", "BC_25", "BC_3", "BC_5")
Naive_wo = c("BC_1", "BC_10", "BC_11", "BC_12", 'BC_13',
             'BC_18', 'BC_4', 'BC_6', 'BC_7', 'BC_9')
Recurrence = c("BC_14", "BC_15", "BC_16", "BC_17", "BC_2", "BC_20",
               "BC_21","BC_22","BC_23","BC_24","BC_26","BC_27","BC_8")

summary_df <- sig_df %>%
  group_by(sample) %>%
  summarise(across(ends_with("_score"), \(x) mean(x, na.rm = TRUE))) %>%
  dplyr::mutate(group = case_when(
    sample %in% Naive_w ~ "Naive_w",
    sample %in% Naive_wo ~ "Naive_wo",
    sample %in% Recurrence ~ "Recurrence",
    TRUE ~ "Unknown"
  ))

# Remove high basal type samples
summary_df = dplyr::filter(summary_df, !sample %in% c('BC_23','BC_2','BC_10','BC_7'))

# Boxplot
plot_df <- summary_df %>%
  pivot_longer(cols = ends_with("_score"), 
               names_to = "signature", 
               values_to = "mean_score")

my_comparisons <- list( 
  c("Naive_wo", "Naive_w"), 
  c("Naive_wo", "Recurrence") 
)

plot_df$group = factor(plot_df$group, levels = c('Naive_wo','Naive_w','Recurrence'))

p = dplyr::filter(plot_df, signature == 'pEMT_pancancer_score') %>%
  ggplot(aes(x = group, y = mean_score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5, shape = 16) +
  facet_wrap(~signature, scales = "free_y") +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    method.args = list(alternative = "less"),
    size = 2,
    label = "p.format", 
    hide.ns = FALSE 
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Pairwise Wilcoxon Test (BH-adjusted)",
    y = "Mean Score per Sample",
    x = ""
  ) +
   theme_classic(base_size = 8) +
   theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

pdf("Reproducibility/Results/Plots/ECOTYPE/Figure7G.pdf", w=5, h=5)
 plot(p)
dev.off()

#========================================================
## IFN

plot_df_2groups <- dplyr::filter(plot_df, group %in% c("Naive_w", "Naive_wo"))

my_comparisons <- list( c("Naive_w", "Naive_wo") )

p = dplyr::filter(plot_df_2groups, signature == 'Interferon_pancancer_score') %>%
  ggplot(aes(x = group, y = mean_score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2, shape = 16) +
  facet_wrap(~signature, scales = "free_y") +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    size = 2,
    label = "p.format") +
  scale_fill_manual(values = c("Naive_w" = "#fc8d62", "Naive_wo" = "#66c2a5")) +
  labs(
    title = "Comparison: Naive_w vs Naive_wo",
    y = "Mean Score per Sample",
    x = ""
  ) +
  theme_classic(base_size = 8) +
  theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

pdf("Reproducibility/Results/Plots/ECOTYPE/FigureS13L.pdf", w=3, h=3)
 plot(p)
dev.off()

#########################################################
## Fig.7H  CD4_T_CD26Hi proportion boxplot 
#########################################################
# External cohort
meta_npj = fread_n("Reproducibility/Results/scANVI/BCG/annotation_transfer_metadata.txt") %>%
            dplyr::select(sample, predicted_celltype)
tmp_tab <- table(meta_npj$sample, meta_npj$predicted_celltype)
cell_counts <- as.data.frame.matrix(tmp_tab)

df_prop = data.frame(100*cell_counts/rowSums(cell_counts)) %>% as.data.frame() %>% t()  # celltype*pts
df_prop = t(df_prop) %>% as.data.frame()
df_prop <- df_prop %>%
  dplyr::mutate(
    iCAF = iCAF_CD321 + iCAF_SLC14A1 + proCAF,
    NK_CD56_CD49a_Hi = NK_CD56_CD49a_Hi_CD103_Hi + NK_CD56_CD49a_Hi_CD103_Lo,
    CD8_Tex = CD8_Tex_1 + CD8_Tex_2
  ) %>%
  dplyr::select(-proCAF, -iCAF_CD321, -iCAF_SLC14A1, -CD8_Tex_1, -CD8_Tex_2,
         -NK_CD56_CD49a_Hi_CD103_Hi, -NK_CD56_CD49a_Hi_CD103_Lo)

cd4_subsets <- c(
  "CD4_CTL", "CD4_T_CD26", "CD4_T_proliferative", "CD4_Tcm",
  "CD4_Tfh.like", "CD4_Th17", "CD4_Tn", "CD4_Tph.like", "CD4_Tsen",'Treg_naive','Treg_effector',
  "CD8_T_proliferative", "CD8_Tcm", "CD8_Tem", "CD8_Temra", "CD8_Tn", "CD8_Trm","CD8_Tex" 
)

df_prop$CD26_frac <- df_prop$CD4_T_CD26 / rowSums(df_prop[, cd4_subsets], na.rm = TRUE)
Recurrence <- c("BC_14", "BC_15", "BC_16", "BC_17", "BC_2", "BC_20", "BC_21","BC_22","BC_23","BC_24","BC_26","BC_27","BC_8")
bcg_prop_npj = df_prop[Recurrence,]
bcg_prop_npj$rec = 'rec'
bcg_prop_npj2 = dplyr::filter(bcg_prop_npj, Epithelial<94)

#################################################
# DOGMA cohort
metadata = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_metadata.txt")

df_immune = dplyr::filter(metadata, lineage %in% c('B','CD4_T','CD8_T_NK_ILC','Myeloid'))
res_immune = as.data.frame.matrix(table(df_immune$sample, df_immune$celltype))
df_prop_DOGMA <- res_immune %>%
  { 100 * . / rowSums(.) } %>%
  as.data.frame()

cd4_subsets2 <- c(
  "CD4_CTL", "CD4_T_CD26", "CD4_T_proliferative", "CD4_Tcm",
  "CD4_Tfh-like", "CD4_Th17", "CD4_Tn", "CD4_Tph-like", "CD4_Tsen", 'Treg_naive','Treg_effector',
  "CD8_T_proliferative", "CD8_Tcm", "CD8_Tem", "CD8_Temra", "CD8_Tn", "CD8_Trm", "CD8_Tex_1", "CD8_Tex_2" 
)
df_prop_DOGMA$CD26_frac <- df_prop_DOGMA$CD4_T_CD26 / rowSums(df_prop_DOGMA[, cd4_subsets2], na.rm = TRUE)

post = c('BC_022','BC_027','BC_032',"BC_039","BC_040",'BC_043',"BC_044","BC_047","BC_048")
rec_dogma = c('rec', 'no_rec', 'rec', 'no_rec','no_rec', 'no_rec', 'no_rec', 'no_rec', 'no_rec')
bcg_prop_DOGMA = df_prop_DOGMA[post,]
bcg_prop_DOGMA$rec = rec_dogma

#################################################
combined_rec_data <- bind_rows(
  bcg_prop_npj2 %>% dplyr::select(CD26_frac, rec) %>% dplyr::mutate(Cohort = "npj"),
  bcg_prop_DOGMA %>% dplyr::select(CD26_frac, rec) %>% dplyr::mutate(Cohort = "DOGMA")
)

combined_rec_data$rec <- factor(combined_rec_data$rec, levels = c("no_rec", "rec"))
my_comparisons <- list(c("no_rec", "rec"))

npg_cols <- pal_npg()(2)
rev_npg_cols <- rev(npg_cols)

p <- ggplot(combined_rec_data, aes(x = rec, y = CD26_frac, fill = rec)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5, shape = 16) + 
  # Wilcoxon検定の追加
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    size = 3.5,
    label = "p.format" 
  ) +
  scale_fill_manual(values = rev_npg_cols) +
  labs(
    title = "CD26_frac by Recurrence Status",
    subtitle = "Combined Cohort (npj + DOGMA)",
    x = "Recurrence",
    y = "CD26+ CD4 T cell Fraction"
  ) +
  theme_classic(base_size = 8) +
  theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

pdf("Reproducibility/Results/Plots/ECOTYPE/Figure7H.pdf", w=3, h=3)
  print(p)
dev.off()

