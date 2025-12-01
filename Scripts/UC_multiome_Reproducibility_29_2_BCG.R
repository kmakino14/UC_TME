#####################################################
# DOGMA analysis - BCG -
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

suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsci))
suppressMessages(library(colorspace))
suppressMessages(library(BuenColors))

suppressMessages(library(rstatix))
suppressMessages(library(dplyr))
suppressMessages(library(ggpubr))

set.seed(1234)

#******************************************
data_dir = "Reproducibility/Data"
lineage = 'BCG'

DOGMA = file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")) %>% readRDS()

DOGMA$prepos = fct_collapse(DOGMA$sample,
                            pre = c("BC_011","BC_016","BC_023","BC_033","BC_037"),
                            post = c('BC_027','BC_032',"BC_039","BC_040",'BC_043',"BC_044","BC_047","BC_048")
                            ) %>% factor(., levels = c("pre","post"))

#################################################
# Fig.5B  Bar plot of differential genes/peaks
#################################################
# RNA
celltypes = c("CD4_Tconv","Treg","CD8_T","NK_ILC","B","Mono_Mac","DC")

deg_counts <- list()
for (tmp_celltype in celltypes) {
  tmp_path <- paste0("Reproducibility/Results/Differential/BCG/UC_DOGMA_DEG_MAST_BCG_", tmp_celltype, ".txt")

  # load
  de_df <- fread_n(tmp_path)
  de_df$DEG <- "nonDEG"
  de_df$DEG[(de_df$avg_log2FC > 0.5 | de_df$avg_log2FC < -0.5) & de_df$p_val_adj < 0.05] <- "DEG"

  # count
  counts <- table(de_df$DEG)

  # ensure both categories exist
  counts_full <- c(DEG = 0, nonDEG = 0)
  counts_full[names(counts)] <- counts

  # store
  deg_counts[[tmp_celltype]] <- counts_full
}

# combine into dataframe
deg_counts <- as.data.frame(deg_counts) %>% t() %>% as.data.frame()
deg_counts$DEG    <- as.integer(deg_counts$DEG)
deg_counts$nonDEG <- as.integer(deg_counts$nonDEG)
deg_counts = dplyr::arrange(deg_counts, DEG)

paste0("Reproducibility/Results/Plots/BCG/Figure5B_barplot_DEG_counts.pdf")%>% pdf(., w = 3.5, h = 5)
 barplot(deg_counts$DEG,
        names.arg = rownames(deg_counts),
        las = 2, col = "#F8766D",xlim = c(0, 1000), horiz = TRUE, 
        ylab = "Number of DEGs", xlab = "Cell type")
dev.off()


#################################################
# ATAC
celltypes = c("CD4_Tconv","Treg","CD8_T","NK_ILC","B","Mono_Mac","DC")

da_counts <- list()
for (tmp_celltype in celltypes) {
    tmp_path = paste0('Reproducibility/Results/SnapATAC2/BCG/snapatac2_BCG_diff_peaks_',tmp_celltype,'.txt')

    # load
    de_df = fread_n(tmp_path)
    de_df$DA = "nonDA"
    de_df$DA[(de_df$log2FC > 0.25 | de_df$log2FC < -0.25) & de_df$p_value_adj < 0.1] <- "DA"

    # count
    counts <- table(de_df$DA)
  
    # ensure both categories exist
    counts_full <- c(DA = 0, nonDA = 0)
    counts_full[names(counts)] <- counts
  
    # store
    da_counts[[tmp_celltype]] <- counts_full
}

# combine into dataframe
da_counts <- as.data.frame(da_counts) %>% t() %>% as.data.frame()
da_counts$DA    <- as.integer(da_counts$DA)
da_counts$nonDA <- as.integer(da_counts$nonDA)
da_counts = dplyr::arrange(da_counts, DA)

paste0("Reproducibility/Results/Plots/BCG/Figure5B_barplot_DA_counts.pdf")%>% pdf(., w = 3.5, h = 5)
 barplot(da_counts$DA,
        names.arg = rownames(deg_counts),
        las = 2, col = "#7cae00",xlim = c(0, 25000), horiz = TRUE, 
        ylab = "Number of DAs", xlab = "Cell type")
dev.off()


#################################################
# Fig.5C  Differential volcano plot
#################################################
path = "Reproducibility/Results/MultiNicheNet/"
multinichenet_output = paste0(path, "multinichenet_output.rds") %>% readRDS()

contrasts_oi = c("'post-pre','pre-post'")
contrast_tbl = tibble(contrast =
                        c("post-pre", "pre-post"),
                      group = c("post", "pre"))

deg_CD4_Tconv = as.data.frame(multinichenet_output$celltype_de) %>%
                dplyr::filter(., cluster_id %in% c('CD4_Tconv')) %>%
                dplyr::filter(., contrast %in% c('post-pre')) %>%
                column_to_rownames('gene') %>%
                dplyr::select(., c('logFC','p_adj'))
colnames(deg_CD4_Tconv) = c('avg_log2FC','p_val_adj')

label_up = c('CXCR6','CCR5','CCR6','CCL20','RORC','DPP4',
             'CD40LG', 'GZMA')
label_down = c('IL23A','TCF7','CCR7','BACH2')


paste0("Reproducibility/Results/Plots/BCG/Figure5C_volcano_CD4_Tconv.pdf") %>% pdf(.,width=3,height=3)
  Volcano_DEG_label_v2(data = deg_CD4_Tconv, celltype = "CD4_Tconv", label_up = label_up, label_down = label_down, 
                       log2FC = 0.5, p_val_thresh = 0.05, ylim = 3.5, thresh = 10^(-3.5))
dev.off()

###############################
deg_Mono_Mac = as.data.frame(multinichenet_output$celltype_de) %>%
                dplyr::filter(., cluster_id %in% c('Mono_Mac')) %>%
                dplyr::filter(., contrast %in% c('post-pre')) %>%
                column_to_rownames('gene') %>%
                dplyr::select(., c('logFC','p_adj')) %>%
                dplyr::arrange(., desc(logFC))
colnames(deg_Mono_Mac) = c('avg_log2FC','p_val_adj')

label_up = c('VCAM1','CXCL9','CXCL10','CXCL12','STAT1')
label_down = c('ITGA5','GPR183','CD109')

paste0("Reproducibility/Results/Plots/BCG/Figure5C_volcano_Mono_Mac.pdf") %>% pdf(.,width=3,height=3)
  Volcano_DEG_label_v2(data = deg_Mono_Mac, celltype = "Mono_Mac", label_up = label_up, label_down = label_down, 
                       log2FC = 0.5, p_val_thresh = 0.05, ylim = 3, thresh = 1e-3)
dev.off()

################################
# Extract gene names
genes_CD4_Tconv <- deg_CD4_Tconv %>%
  dplyr::filter(avg_log2FC > 0.5 | p_val_adj < 0.05) %>%
  rownames()

genes_Mono_Mac <- deg_Mono_Mac %>%
  dplyr::filter(avg_log2FC > 0.5 | p_val_adj < 0.05) %>%
  rownames()

# Create data frame for Metascape format
metascape_df <- data.frame(
  Name  = c("CD4_Tconv", "Mono_Mac"),
  Genes = c(
    paste(genes_CD4_Tconv, collapse = ","),
    paste(genes_Mono_Mac, collapse = ",")
  )
)

# Write to tab-separated .txt file
write.table(
  metascape_df,
  file = "Reproducibility/Results/MultiNicheNet/post_BCG_gene_list_for_metascape.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#################################################
# Fig.5C/5D GO plot
#################################################
# Metascape
GO_df  = fread("Reproducibility/Results/MultiNicheNet/GO_AllLists.csv") %>% as.data.frame() %>%
         dplyr::select("GO","Description","Log(q-value)","GeneList")
colnames(GO_df)[3] = "FDR" 
GO_df$FDR = -(GO_df$FDR)

# CD4 Tconv
GO_term_use = c("GO:0001819","GO:0050900","GO:0042110","GO:0022409","GO:0030217","GO:0006954","GO:0032729")
GO_result = dplyr::filter(GO_df, GeneList %in% c("CD4_Tconv")) %>% 
            dplyr::filter(GO %in% GO_term_use) %>% 
            dplyr::distinct(., GO,.keep_all = TRUE) %>% 
            arrange(., desc(FDR))
GO_result$Description = factor(GO_result$Description, levels = rev(GO_result$Description))

p1 = ggplot(GO_result, aes(x = Description, y= FDR, label = Description)) + 
     geom_bar(stat = "identity",fill = "#0081c9") +
     coord_flip() + labs(x = "", y = "-log(FDR)") +
     theme_bw()+
     theme(panel.grid.minor = element_blank())

paste0("Reproducibility/Results/Plots/BCG/Figure5C_barplot_CD4_Tconv.pdf") %>% pdf(.,h=5,w=15)
 plot(p1)
dev.off()  

# Mono_Mac
GO_term_use = c("GO:0050778","GO:0006954","GO:0001819","GO:0009617","GO:0071345","GO:0032946","GO:0034341")
GO_result = dplyr::filter(GO_df, GeneList %in% c("Mono_Mac")) %>% 
            dplyr::filter(GO %in% GO_term_use) %>% 
            dplyr::distinct(., GO,.keep_all = TRUE) %>% 
            arrange(., desc(FDR))
GO_result$Description = factor(GO_result$Description, levels = rev(GO_result$Description))

p2 = ggplot(GO_result, aes(x = Description, y= FDR, label = Description)) + 
     geom_bar(stat = "identity",fill = "#ff5a00") +
     coord_flip() + labs(x = "", y = "-log(FDR)") +
     theme_bw()+
     theme(panel.grid.minor = element_blank())

paste0("Reproducibility/Results/Plots/BCG/Figure5D_barplot_Mono_Mac.pdf") %>% pdf(.,h=5,w=15)
 plot(p2)
dev.off()


#################################################
# Fig.5E  Violin plot of signature scores
#################################################

DOGMA$coarse_celltype_w_BCG = paste0(DOGMA$coarse_celltype, "_", DOGMA$prepos)

Idents(DOGMA) = 'coarse_celltype_w_BCG'
DOGMA_tmp = subset(DOGMA, ident = c('Mono_Mac_pre', 'Mono_Mac_post'))
sig_df = fread_n('Reproducibility/Results/VISIONR/UC_DOGMA_BCG_signature_score_literature.txt')
DOGMA_tmp[["signature"]] = CreateAssayObject(data = sig_df %>% t() %>% as.matrix())

DefaultAssay(DOGMA_tmp) <- "signature"
sig_mat <- DOGMA_tmp[["signature"]]@data
genes   <- rownames(sig_mat)
groups  <- DOGMA_tmp@meta.data$coarse_celltype_w_BCG

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
  dplyr::filter(Group %in% c("Mono_Mac_pre", "Mono_Mac_post")) %>%
  mutate(
    Group = factor(Group, levels = c("Mono_Mac_pre", "Mono_Mac_post")),
    Gene  = factor(Gene, levels = genes)  # facet order = rownames(sig assay)
  )

# Colors
colors <- c("Mono_Mac_pre" = "#2461a1", "Mono_Mac_post" = "#a91e2c")

# Pairwise comparison
my_comparisons <- list(c("Mono_Mac_pre", "Mono_Mac_post"))

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
  facet_wrap(~ Gene, scales = "free_y") + 
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "fdr",
    bracket.nudge.y = 0.1,
    hide.ns = TRUE,
    tip.length = 0
  )

pdf_3("Reproducibility/Results/Plots/BCG/Figure5E.pdf", h = 4, w = 9)
 plot(p)
dev.off()


#################################################
# Fig.5F  Milo beeswarm plot
#################################################
df = fread_n(paste0("Reproducibility/Results/Milo/output/Milo_BCG_design_patient_BCG_contrasts_post_vs_pre_result_TAM_merge_ver.txt"))

p = df %>%
    dplyr::filter(nhood_annotation_frac > 0.5) %>%      ## 0.35にするとCD8_Trmが見えるが…
    mutate(signif = ifelse(SpatialFDR < 0.1, logFC, 0)) %>%
    group_by(nhood_annotation) %>%
    mutate(mean_lfc_val = ifelse(!is.na(signif), logFC, 0)) %>%
    mutate(mean_lfc = mean(mean_lfc_val)) %>%
    ungroup() %>%
    dplyr::arrange(mean_lfc) %>%
    mutate(nhood_annotation = factor(nhood_annotation, levels = unique(nhood_annotation))) %>%
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
                          limits = c(-6, 6)) +
    theme_bw(base_size = 18) +
    geom_hline(yintercept = 0, linetype = 2) +
    xlab('') + 
    ylim(-6, 6)

paste0("Reproducibility/Results/Plots/BCG/Figure5F.pdf") %>% pdf(.,width=11,height=12)
 plot(p)
dev.off()


#################################################
# Fig.S10B  DPP4/IFNG/CD26 expression
#################################################
Idents(DOGMA) = "sample"
BCG_paired = subset(DOGMA, ident = c("BC_011","BC_039","BC_023","BC_044","BC_033","BC_048","BC_037","BC_047")) 
BCG_paired$timepoint = fct_collapse(BCG_paired$sample, 
                             pre  = c("BC_011","BC_023","BC_033","BC_037"),
                             post = c("BC_039","BC_044","BC_048","BC_047")
                             )
BCG_paired$patient = fct_collapse(BCG_paired$sample, 
                           P1 = c("BC_011","BC_039"),
                           P2 = c("BC_023","BC_044"),
                           P3 = c("BC_033","BC_048"),
                           P4 = c("BC_037","BC_047")
                           )
BCG_paired$sample_id = paste0(BCG_paired$patient,"_",BCG_paired$timepoint)


Idents(BCG_paired) = 'coarse_celltype'
BCG_paired = subset(BCG_paired, ident = c('CD4_Tconv')) 

RNA_pseudobulk = AverageExpression(BCG_paired, assay = "RNA", slot = "data", group.by = "sample_id", return.seurat = TRUE) 
ADT_pseudobulk = AverageExpression(BCG_paired, assay = "ADT", slot = "data", group.by = "sample_id", return.seurat = TRUE) 

RNA_df_long <- as.data.frame(RNA_pseudobulk[["RNA"]]$data) %>% rownames_to_column('gene') %>% 
  dplyr::filter(gene %in% c('DPP4','IFNG','GZMA','GZMH')) %>% 
  pivot_longer(
    cols = -c('gene'),
    names_to = "sample_id",
    values_to = "value"
    ) %>%
    mutate(
      group = if_else(str_detect(sample_id, "pre"), "pre", "post"),
      pair = take_factor(sample_id,1,"_")
    ) %>% mutate(group = factor(group, levels = c("pre", "post"))) %>%
    arrange(gene, pair, group)

ADT_df_long <- as.data.frame(ADT_pseudobulk[["ADT"]]$data) %>% rownames_to_column('gene') %>% 
  dplyr::filter(gene %in% c('surface-A0396-CD26')) %>% 
  pivot_longer(
    cols = -c('gene'),
    names_to = "sample_id",
    values_to = "value"
    ) %>%
    mutate(
      group = if_else(str_detect(sample_id, "pre"), "pre", "post"),
      pair = take_factor(sample_id,1,"_")
    ) %>% mutate(group = factor(group, levels = c("pre", "post"))) %>%
    arrange(gene, pair, group)

RNA_df_long$data_type = 'RNA'
ADT_df_long$data_type = 'ADT'
combined_df <- bind_rows(RNA_df_long, ADT_df_long)

p = ggplot(combined_df, aes(x = group, y = value)) +
  geom_line(aes(group = pair), color = "black", size = 1) +
  geom_point(aes(fill = group), shape = 21, size = 4, color = "black", stroke = 0.5) +
  # add pair labels on the "post" side
  geom_text(
    data = combined_df %>% dplyr::filter(group == "post"),
    aes(label = pair),
    hjust = -0.2, vjust = 0.5, size = 3.5
  ) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_classic() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlim(c("pre","post")) +
  scale_fill_manual(
    values = c("pre" = "#245d9880", "post" = "#a1212c80")
  )

paste0("Reproducibility/Results/Plots/BCG/FigureS10B.pdf") %>% pdf(., w=7, h=4)
 plot(p)
dev.off()


#################################################
# Fig.S10C  CD4_T_CD26 volcano plot
#################################################
data = fread_n("scvi/CD4/library_ver0.18.0_w_ADT_select/DEG/UC_DOGMA_reseq_deg_MAST_CD4_T_CD26_by_status.txt")
data$avg_log2FC[data$avg_log2FC > 4] <- 4
data$avg_log2FC[data$avg_log2FC < -4] <- -4
label_up = c('RUNX1','DPP4','HDAC9','STAT1','CD38','TET1')
label_down = c('IL22','CCL2','CCL3','CCL4','IL17A','CXCL8')

paste0("Reproducibility/Results/Plots/BCG/FigureS10B.pdf") %>% pdf(.,width=3,height=3)
  Volcano_DEG_label_v2(data = data, celltype = "status", 
     label_up = label_up, label_down = label_down, log2FC = 0.5, ylim = 100, thresh = 1e-100)
dev.off()

#################################################
# Fig.S10D  MAESTER
#################################################
source("Reproducibility/Scripts/Source/variant_calling.R")

maegatk.rse = readRDS("Reproducibility/Data/MAESTER/maegatk.rse.rds")
af.dm.all = readRDS("Reproducibility/Data/MAESTER/af_dm_all.rds")
metadata_list = readRDS("Reproducibility/Data/MAESTER/metadata_list.rds")

grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
blue_scale <- brewer.pal(n= 9,name="Blues")
colorscale3 <- c(grey_scale[3], grey_scale[3], grey_scale[3], red_scale[2:9])
colorscale4 <- c(grey_scale[3], grey_scale[3], grey_scale[3], blue_scale[2:9])

#########################################################
sample_before = c("BC_011","BC_023","BC_033","BC_037")
sample_after  = c("BC_039","BC_044","BC_048","BC_047")
v_list = c("5220_C>T",'6048_G>A','11969_G>A','1787_G>A')

plots <- list()
for(tmp_num in 1:4){
    tmp_sample_before = sample_before[tmp_num]
    tmp_sample_after  = sample_after[tmp_num]
    v = v_list[tmp_num]

    # Intersect barcodes for maegatk and metadata
    tmp_metadata <- metadata_list[[tmp_num]] %>% as.data.frame() %>% group_by(cell) %>% dplyr::filter(n() == 1) %>% ungroup
    tmp_metadata.tib <- tmp_metadata %>% dplyr::filter(cell %in% colnames(af.dm.all))
    tmp_maegatk <- maegatk.rse[,tmp_metadata.tib$cell]
    af.dm = af.dm.all[,tmp_metadata.tib$cell]

    CellSubsets.ls <- list(unionCells = tmp_metadata.tib$cell)
    lengths(CellSubsets.ls)

    # Add info of variant of interest
    tmp_metadata.tib$cov_voi <- assays(tmp_maegatk)[["coverage"]][as.numeric( cutf(v, d = "_") ),tmp_metadata.tib$cell]

    tmp_metadata.tib$af_voi <- af.dm[v,tmp_metadata.tib$cell]
    tmp_metadata.tib$af_voi[tmp_metadata.tib$af_voi<2] <- 0
    tmp_metadata.tib$af_voi_pre = tmp_metadata.tib$af_voi
    tmp_metadata.tib$af_voi_pre[tmp_metadata.tib$sample==tmp_sample_after] = 0
    tmp_metadata.tib$af_voi_post = tmp_metadata.tib$af_voi
    tmp_metadata.tib$af_voi_post[tmp_metadata.tib$sample==tmp_sample_before] = 0

    tmp_metadata.tib_pre = dplyr::filter(tmp_metadata.tib, !STAGE %in% c("post_BCG"))
    tmp_metadata.tib_post = dplyr::filter(tmp_metadata.tib, STAGE %in% c("post_BCG"))
 
    tmp_post = dplyr::filter(tmp_metadata.tib_post, af_voi_post>1)

    p1 =    tmp_metadata.tib %>% dplyr::arrange(af_voi_pre) %>% #dplyr::filter(cov_voi > 5) %>%
            ggplot(aes(x = UMAP_1, y = UMAP_2, color = af_voi_pre)) + # change to cov_voi to see coverage
            geom_point_rast(size = 1) +
            scale_color_gradientn(colors = colorscale4, limits = c(0,100), n.breaks = 3) +
            theme_classic() +
            theme(aspect.ratio = 1, axis.line = element_blank(), plot.title = element_text(hjust=0.5),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            ggtitle(v)

    p2 =    tmp_metadata.tib %>% dplyr::arrange(af_voi_post) %>% #dplyr::filter(cov_voi > 5) %>%
            ggplot(aes(x = UMAP_1, y = UMAP_2, color = af_voi_post)) + # change to cov_voi to see coverage
            geom_point_rast(size = 1) +
            scale_color_gradientn(colors = colorscale3, limits = c(0,100), n.breaks = 3) +
            theme_classic() +
            theme(aspect.ratio = 1, axis.line = element_blank(), plot.title = element_text(hjust=0.5),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
            ggtitle(v)

    plots[[tmp_num]] <- p1|p2 

}
paste0("Reproducibility/Results/Plots/BCG/FigureS10D.pdf") %>% pdf(,width = 12, height = 6)
 plots
dev.off()


#################################################
# Fig.5J  Scatter plot of TF expression/activity
#################################################
# CD4_Tconv
TF_exp_df = fread_n('Reproducibility/Results/LINGER/BCG/cell_population_exp_TF_activity_zscore_CD4_Tconv_BCG_paired.txt') %>%
            rownames_to_column('Name')

highlight_pos <- c("STAT1", "RORA", "RUNX2", "RUNX1", 'FOSL2')
highlight_neg <- c('BACH2', 'TCF7')

# Modify the ggplot code
p <- ggplot(TF_exp_df, aes(x = EXP, y = TF)) +
    geom_point(data = subset(TF_exp_df, !(Name %in% c(highlight_pos, highlight_neg))),
               aes(color = "normal"), alpha = 0.5, shape = 16, size = 1) +
    geom_point(data = subset(TF_exp_df, Name %in% highlight_pos),
               aes(fill = "highlight_pos"), shape = 21, size = 2) +
    geom_point(data = subset(TF_exp_df, Name %in% highlight_neg),
               aes(fill = "highlight_neg"), shape = 21, , size = 2) +
    geom_text_repel(data = subset(TF_exp_df, Name %in% c(highlight_pos, highlight_neg)),
                    aes(label = Name, color = case_when(
                        Name %in% highlight_pos ~ "highlight_pos",
                        Name %in% highlight_neg ~ "highlight_neg"
                    )), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("highlight_pos" = "#A92524", "highlight_neg" = "#467DAB")) +
    scale_color_manual(values = c("highlight_pos" = "#A92524", "highlight_neg" = "#467DAB", "normal" = "grey")) +
    theme_classic() +
    labs(title = "CD4_Tconv") +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          panel.border = element_rect(color = "black", fill = NA))

#################################################
# Mono_Mac
TF_exp_df = fread_n('Reproducibility/Results/LINGER/BCG/cell_population_exp_TF_activity_zscore_Mono_Mac_BCG_paired.txt') %>%
            rownames_to_column('Name')

highlight_pos <- c("STAT1", "IRF2", "IRF8", "JDP2")
highlight_neg <- c('HIF1A', 'NFAT5', 'KLF4')

# Modify the ggplot code
q <- ggplot(TF_exp_df, aes(x = EXP, y = TF)) +
    geom_point(data = subset(TF_exp_df, !(Name %in% c(highlight_pos, highlight_neg))),
               aes(color = "normal"), alpha = 0.5, shape = 16, size = 1) +
    geom_point(data = subset(TF_exp_df, Name %in% highlight_pos),
               aes(fill = "highlight_pos"), shape = 21, size = 2) +
    geom_point(data = subset(TF_exp_df, Name %in% highlight_neg),
               aes(fill = "highlight_neg"), shape = 21, , size = 2) +
    geom_text_repel(data = subset(TF_exp_df, Name %in% c(highlight_pos, highlight_neg)),
                    aes(label = Name, color = case_when(
                        Name %in% highlight_pos ~ "highlight_pos",
                        Name %in% highlight_neg ~ "highlight_neg"
                    )), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("highlight_pos" = "#A92524", "highlight_neg" = "#467DAB")) +
    scale_color_manual(values = c("highlight_pos" = "#A92524", "highlight_neg" = "#467DAB", "normal" = "grey")) +
    theme_classic() +
    labs(title = "Mono_Mac") +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          panel.border = element_rect(color = "black", fill = NA))

paste0("Reproducibility/Results/Plots/BCG/Figure5J.pdf") %>% pdf(., w=7, h=3.5)
 plot(p|q)
dev.off()

#################################################
# Fig.S10F  Violin plot of TF activity
#################################################
# CD4_T
CD4_TF_df = fread_n("Reproducibility/Results/LINGER/BCG/cell_population_TF_activity_CD4_T_BCG.txt")
cells.use = intersect(colnames(DOGMA), colnames(CD4_TF_df))
DOGMA_tmp = subset(DOGMA, cells = cells.use)
DOGMA_tmp[["TF_activity"]]   = CreateAssayObject(data = CD4_TF_df[,cells.use] %>% as.matrix())

DefaultAssay(DOGMA_tmp) <- 'TF_activity'
DOGMA_tmp = ScaleData(DOGMA_tmp, scale.max = 10)

Idents(DOGMA_tmp) = 'coarse_celltype'
DOGMA_tmp = subset(DOGMA_tmp, ident = c('CD4_Tconv'))
target <- c('STAT1', 'RORA', 'RUNX1')
combined_df <- data.frame()

# Loop through the target genes to prepare the data
for(tmp_gene in target){
  tmp_df <- data.frame(
    Group = DOGMA_tmp@meta.data$prepos,
    Value = DOGMA_tmp[['TF_activity']]$scale.data[tmp_gene, ],
    Gene = tmp_gene
  )
  combined_df <- rbind(combined_df, tmp_df)
}

# Set the group factor levels
combined_df$Group <- factor(combined_df$Group, levels = c('pre','post'))
combined_df$Gene <- factor(combined_df$Gene, levels = target)

# Define custom colors
colors <- c('pre' = "#2461a1", 'post' = '#a91e2c')

# Create the plot
p <- ggplot(combined_df, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with fill by group
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +  # Boxplot overlay
  scale_fill_manual(values = colors) +  # Apply custom colors
  theme_classic() +
  labs(y = "Score") +
  theme(legend.position = "none") +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +  # Facet by gene+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +  # Add margin to the y-axis
  stat_compare_means(
    comparisons = list(c('pre', 'post')),
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "BH",
    bracket.nudge.y = 0.1,
    hide.ns = TRUE,
    tip.length = 0
  )

paste0("Reproducibility/Results/Plots/BCG/FigureS10F_CD4_Tconv.pdf") %>% pdf_3(.,h=3.5,w=7)
 plot(p)
dev.off()

#################################################
# Mono_Mac
Myeloid_TF_df = fread_n("Reproducibility/Results/LINGER/BCG/cell_population_TF_activity_Myeloid_BCG.txt")
cells.use = intersect(colnames(DOGMA), colnames(Myeloid_TF_df))
DOGMA_tmp = subset(DOGMA, cells = cells.use)
DOGMA_tmp[["TF_activity"]]   = CreateAssayObject(data = Myeloid_TF_df[,cells.use] %>% as.matrix())

DefaultAssay(DOGMA_tmp) <- 'TF_activity'
DOGMA_tmp = ScaleData(DOGMA_tmp, scale.max = 10)

Idents(DOGMA_tmp) = 'coarse_celltype'
DOGMA_tmp = subset(DOGMA_tmp, ident = c('Mono_Mac'))
target <- c('STAT1', 'IRF8', 'JDP2')
combined_df <- data.frame()

# Loop through the target genes to prepare the data
for(tmp_gene in target){
  tmp_df <- data.frame(
    Group = DOGMA_tmp@meta.data$prepos,
    Value = DOGMA_tmp[['TF_activity']]$scale.data[tmp_gene, ],
    Gene = tmp_gene
  )
  combined_df <- rbind(combined_df, tmp_df)
}

# Set the group factor levels
combined_df$Group <- factor(combined_df$Group, levels = c('pre','post'))
combined_df$Gene <- factor(combined_df$Gene, levels = target)

# Define custom colors
colors <- c('pre' = "#2461a1", 'post' = '#a91e2c')

# Create the plot
p <- ggplot(combined_df, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with fill by group
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +  # Boxplot overlay
  scale_fill_manual(values = colors) +  # Apply custom colors
  theme_classic() +
  labs(y = "Score") +
  theme(legend.position = "none") +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +  # Facet by gene+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +  # Add margin to the y-axis
  stat_compare_means(
    comparisons = list(c('pre', 'post')),
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "BH",
    bracket.nudge.y = 0.1,
    hide.ns = TRUE,
    tip.length = 0
  )

paste0("Reproducibility/Results/Plots/BCG/FigureS10F_Mono_Mac.pdf") %>% pdf_3(.,h=3.5,w=7)
 plot(p)
dev.off()


#########################################################
## Fig.5K/S10G/S10H  Gene track with co-accessible regions
#########################################################
source("Reproducibility/Scripts/Source/Seurat_source.R")

DefaultAssay(DOGMA) = 'ATAC'
Idents(DOGMA) = "sample"
BCG_paired = subset(DOGMA, ident = c("BC_011","BC_039","BC_023","BC_044","BC_033","BC_048","BC_037","BC_047")) 

#########################################################
# CD4_Tconv
Idents(BCG_paired) = "coarse_celltype"
CD4_Tconv = subset(BCG_paired, ident = c('CD4_Tconv'))

region_list <- list(
  CCL20 = c("chr2-227761000-227850000"),
  DPP4 = c("chr2-161948000-162080000")
)

motif_list = list(
  CCL20 = c('RORA','FOSL2','STAT1'),
  DPP4  = c('STAT1','FOSL2')
)

cols <- c('#245D98', '#A1212C')

Idents(CD4_Tconv) = "prepos"
Plot_genetrack_BCG(Obj=CD4_Tconv,lineage='CD4_T',group_list=c('pre','post'),gene='CCL20',
               region_list=region_list, motif_list=motif_list,cutoff=0.2,cols=cols,
               path='Reproducibility/Results/Plots/BCG/Figure5K_CCL20.pdf')

Plot_genetrack_BCG(Obj=CD4_Tconv,lineage='CD4_T',group_list=c('pre','post'),gene='DPP4',
               region_list=region_list, motif_list=motif_list,cutoff=0.1,cols=cols,
               path='Reproducibility/Results/Plots/BCG/FigureS10G_DPP4.pdf')

#########################################################
# Mono_Mac
Idents(BCG_paired) = "coarse_celltype"
Mono_Mac = subset(BCG_paired, ident = c('Mono_Mac'))

region_list <- list(
  CXCL9 = c("chr4-75990000-76050000"),
  CXCL10 = c("chr4-75990000-76050000")
)

motif_list = list(
  CXCL9 = c('STAT1'),
  CXCL10 = c('STAT1')
)

cols <- c('#245D98', '#A1212C')

Idents(Mono_Mac) = "prepos"
Plot_genetrack_BCG(Obj=Mono_Mac,lineage='Myeloid',group_list=c('pre','post'),gene='CXCL9',
               region_list=region_list, motif_list=motif_list,cutoff=0.1,cols=cols,
               path='Reproducibility/Results/Plots/BCG/FigureS10H_CXCL9.pdf')

Plot_genetrack_BCG(Obj=Mono_Mac,lineage='Myeloid',group_list=c('pre','post'),gene='CXCL10',
               region_list=region_list, motif_list=motif_list,cutoff=0.1,cols=cols,
               path='Reproducibility/Results/Plots/BCG/FigureS10H_CXCL10.pdf')


