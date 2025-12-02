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
suppressMessages(library(corrplot))

suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsci))
suppressMessages(library(colorspace))
suppressMessages(library(scales))

suppressMessages(library(dplyr))
suppressMessages(library(ggpubr))

set.seed(1234)

Blue_Red = rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7", "#FFFFFF", "#D1E5F0",
                 "#92C5DE", "#4393C3", "#2166AC", "#053061"))
Orange_Purple = rev(brewer.pal(11, "PuOr"))

#################################################
# Fig.6A  ECOTYPE DOGMA-seq data
#################################################
metadata = fread_n("Reproducibility/Data/UC_DOGMA_metadata.txt")
FACS_prop = fread_n("Reproducibility/Data/UC_DOGMA_FACS_data.txt")

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
                save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure6A_celltype_module.pdf")

#################################################################
df_cor = cor(df_norm2 %>% as.matrix() %>% t() ,method='spearman')

row.order = c('Epithelial', rownames(df_norm2[pheatmap_result$tree_row$order,]))
clustered_data <- df_norm[row.order,]

col_anno_df = fread_n("Reproducibility/Data/TableS1_Clinical_metadata_DOGMA.txt") %>% as.data.frame() %>%
              dplyr::select(., c("status","main_histology"))
col_anno_df$status = factor(col_anno_df$status, levels = c("Early_L","Early_H","Early_UTUC","Advanced_BC","Advanced_UTUC","post_BCG"))
col_anno_df$main_histology = factor(col_anno_df$main_histology,levels = c("UC","Small_cell","Inverted_papilloma","No_malignancy"))
row_anno_df = data.frame(ECOTYPE = rep(c('Tumor','EC1','EC2','EC3','EC4','EC5','EC6','EC7','EC8','EC9', 'EC10'), 
                         times = c(1,4,8,3,7,4,2,7,2,3,2)))
rownames(row_anno_df) = rownames(clustered_data)

# color 
pal1 <- ggsci::scale_color_d3('category10')
pal2 <- ggsci::scale_color_tron()

col.stage <- pal1$palette(6)
col.histo <- pal2$palette(4)
col.eco <- hue_pal()(12)

names(col.stage) <- c("Early_L","Early_H","Early_UTUC","Advanced_BC","Advanced_UTUC","post_BCG")
names(col.histo) <- c("UC","Small_cell","Inverted_papilloma","No_malignancy")
names(col.eco) <- c('Tumor','EC1','EC2','EC3','EC4','EC5','EC6',"EC7","EC8",'EC9','EC10','EC11')

anno_color <- list(
  status = col.stage, 
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
                                    filename=paste0("Reproducibility/Results/Plots/ECOTYPE/Figure6A_celltype_abundance_concat.pdf"), 
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
                   filename=paste0("Reproducibility/Results/Plots/ECOTYPE/Figure6A_celltype_abundance_concat_annotate.pdf"), 
                   border_color = "black" ,na_col = "gray90",
                   main = paste0("ranknorm_spearman_clustered_ecotype"))


##################################################
# Fig.S11A  ECOTYPE TNK external validation
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
                save_path ="Reproducibility/Results/Plots/ECOTYPE/FigureS11A_external.pdf")

########################
## In house dataset
########################
metadata = fread_n("Reproducibility/Data/UC_DOGMA_metadata.txt")
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
                save_path ="Reproducibility/Results/Plots/ECOTYPE/FigureS11A_in_house.pdf")


#################################################
# Fig.6B/Fig.S11B  Signature score boxplot
#################################################
# Load 
sig_df_1 = fread_n("Reproducibility/Results/VISIONR/UC_DOGMA_Malignant_signature_score_literature.txt")
sig_df_2 = fread_n("Reproducibility/Results/VISIONR/UC_DOGMA_Malignant_signature_score_hotspot.txt")
sig_df = cbind(sig_df_1,sig_df_2)
metadata = fread_n("Reproducibility/Data/UC_DOGMA_metadata.txt") %>% dplyr::filter(., celltype %in% c("LUM",'NRP','SQM','MES','NEC'))
metadata = metadata[rownames(sig_df),]

sig_df_w_meta = cbind(sig_df, metadata) %>% dplyr::select(.,c(colnames(sig_df), 'sample', 'STAGE')) 
sig_df_avg <- sig_df_w_meta %>%
  group_by(sample) %>%
  summarise(across(-STAGE, mean, na.rm = TRUE)) %>% as.data.frame()

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
                  tmp_TME = "TME_1", save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure6B_TME_1_vs_others.pdf")
plot_tme_boxplots(df = sig_df_norm, alt = NULL, tmp_feature = "Cytokines_and_inflammatory_response", 
                  tmp_TME = "TME_3", save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure6B_TME_3_vs_others.pdf")
plot_tme_boxplots(df = sig_df_norm, alt = NULL, tmp_feature = "Regulation_of_autophagy", 
                  tmp_TME = "TME_6", save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure6B_TME_6_vs_others.pdf")

plot_tme_boxplots_3group(df = sig_df_norm, alt = 'less', tmp_feature = "Luminal", 
                         save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure6B_LUM.pdf")
plot_tme_boxplots_3group(df = sig_df_norm, alt = 'greater', tmp_feature = "Partial_epithelial-mesenchymal_transition", 
                         save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure6B_pEMT.pdf")


#################################################
# Fig.6C  TF activity boxplot
#################################################
TF_activity = fread_n("LINGER/Primary/Malignant/output/cell_population_TF_activity_zscore.txt")

metadata = fread_n("Reproducibility/Data/UC_DOGMA_metadata.txt") %>% dplyr::filter(., celltype %in% c("LUM",'NRP','SQM','MES','NEC'))
metadata = metadata[rownames(TF_activity),]

TF_df_w_meta = cbind(TF_activity, metadata) %>% 
               dplyr::select(.,c(colnames(TF_activity), 'sample', 'STAGE')) 
TF_df_avg <- TF_df_w_meta %>%
  group_by(sample) %>%
  summarise(across(-STAGE, mean, na.rm = TRUE)) %>% as.data.frame()

TF_df_avg=TF_df_avg %>% column_to_rownames('sample')
TF_df_norm <- as.data.frame(apply(TF_df_avg, 2, function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}))

TF_df_norm = TF_df_norm[c(TME_1,TME_3,TME_4,TME_5,TME_6), ]

plot_tme_boxplots(df = TF_df_norm, alt = NULL, tmp_feature = "ATF4", 
                  tmp_TME = "TME_1", save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure6C_ATF4.pdf")
plot_tme_boxplots(df = TF_df_norm, alt = NULL, tmp_feature = "NR2F2", 
                  tmp_TME = "TME_6", save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure6C_NR2F2.pdf")
plot_tme_boxplots_3group(df = TF_df_norm, alt = 'less', tmp_feature = "GATA3", 
                         save_path = "Reproducibility/Results/Plots/ECOTYPE/Figure6C_GATA3.pdf")

#########################################################
## Fig.6D  Spatial celltype enrichment
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

paste0("Reproducibility/Results/Plots/Slide-seqV2/Figure6D_enrichment.pdf") %>% pdf(w=8, h=10)
 p
dev.off()

#################################################
# Fig.S11F  ECOTYPE Slide-seqV2
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
                   filename=paste0("Reproducibility/Results/Plots/Slide-seqV2/FigureS11F.pdf"), 
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
