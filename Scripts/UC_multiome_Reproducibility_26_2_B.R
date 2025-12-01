#####################################################
# DOGMA analysis - B -
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
lineage = 'B'

DOGMA = file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")) %>% readRDS()

#################################################
# Fig.S6B  Milo beeswarm plot
#################################################
res1 = fread_n("Reproducibility/Results/Milo/output/Milo_B_design_Organ_STAGE_contrasts_Early_H_vs_Early_L_result.txt")
res2 = fread_n("Reproducibility/Results/Milo/output/Milo_B_design_Organ_STAGE_contrasts_Advanced_vs_Early_L_result.txt")

p1 = res1 %>%
     dplyr::filter(nhood_annotation_frac > 0.5) %>%
     mutate(signif=ifelse(SpatialFDR < 0.1, logFC, 0)) %>%
     group_by(nhood_annotation)%>%
     mutate(mean_lfc_val = ifelse(!is.na(signif), logFC, 0)) %>%
     mutate(mean_lfc = mean(mean_lfc_val)) %>%
     ungroup() %>%
     dplyr::arrange(mean_lfc) %>%
     mutate(nhood_annotation=factor(nhood_annotation, 
            levels=c('Plasma',"GC_B",'Atypical_B','B_memory','B_naive'))) %>%
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
            levels=c('Plasma',"GC_B",'Atypical_B','B_memory','B_naive'))) %>%
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

paste0("Reproducibility/Results/Plots/B/FigureS6B.pdf") %>% pdf_3(.,w = 14, h = 4)
 p1+p2
dev.off()


#################################################
# Fig.S6C  Proportion box plot
#################################################
ct_keep <- c("B_naive","B_memory","Atypical_B","GC_B","Plasma")

metadata_IM <- fread_n("Reproducibility/Data/UC_DOGMA_metadata.txt") %>% 
  mutate(celltype = if_else(celltype == "Normal", "Normal", celltype)) %>%
  dplyr::filter(!lineage %in% c("MSC","Endothelial","Epithelial")) %>%
  dplyr::select(celltype, sample, STAGE)

sample_stage <- metadata_IM %>% distinct(sample, STAGE)
cell_frac <- metadata_IM %>%
  dplyr::count(sample, celltype, name = "n") %>%                # counts per (sample, celltype)
  left_join(sample_stage, by = "sample") %>%             # attach STAGE
  group_by(sample, STAGE) %>%
  complete(celltype = ct_keep, fill = list(n = 0)) %>%   # force all 5 types, fill zeros
  mutate(prop = 100 * n / sum(n)) %>%
  ungroup() %>%
  # keep only your 5 cell types (now guaranteed to exist, even if 0)
  dplyr::filter(celltype %in% ct_keep) %>%
  mutate(
    celltype = factor(celltype, levels = ct_keep),
    STAGE    = factor(STAGE, levels = c("Early_L","Early_H","Advanced","post_BCG"))
  )

stage_comparisons <- combn(levels(cell_frac$STAGE), 2, simplify = FALSE)
palette <- c("#00AFBB", "#E7B800", "#FC4E07",'#CC61B0')

# Create boxplot grouped by STAGE_v2 for each variable, with pairwise comparisons
p = ggplot(cell_frac, aes(x = STAGE, y = prop, fill = STAGE)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6, shape = 16) +
  scale_fill_manual(values = palette) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 8) +
  labs(x = "Status", y = "Propotion") +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  geom_pwc(
    aes(group = STAGE),
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "fdr",
    bracket.nudge.y = 0.1,
    hide.ns = TRUE,
    tip.length = 0
  )

paste0("Reproducibility/Results/Plots/B/FigureS6C.pdf") %>% pdf_3(.,w = 14, h = 4)
 p
dev.off()

#################################################
# Fig.S6E  eRegulon heatmap
#################################################

dotplot_df = fread_n("Reproducibility/Results/scenicplus/B/eRegulon_heatmap_dotplot_extended_plotting_df.txt") %>% 
             dplyr::filter(., repressor_activator %in% c('activator'))
colnames(dotplot_df) = c('celltype', 'extended_gene_based_AUC', "eRegulon_name", "extended_region_based_AUC", "repressor_activator")
group = c('B_naive', 'B_memory', 'Atypical_B', 'GC_B', 'Plasma')
dotplot_df$celltype = factor(dotplot_df$celltype, levels = rev(group))

tf_list <- c(
  "KLF2_extended_+/+",
  "ELF2_extended_+/+",
  "PAX5_extended_+/+",
  "SPIB_extended_+/+",
  
  "NFATC2_extended_+/+",
  "RELB_extended_+/+",
  
  "FLI1_extended_+/+",
  "NFAT5_extended_+/+",
  "NFKB1_extended_+/+",
  
  "E2F1_extended_+/+",
  "E2F2_extended_+/+",
  "ETS1_extended_+/+",
  "FOXO1_extended_+/+",
  "IRF8_extended_+/+",
  "YY1_extended_+/+",
  
  "ATF3_extended_+/+",
  "BACH1_extended_+/+",
  "FOS_extended_+/+",
  "IRF4_extended_+/+",
  "JUN_extended_+/+"
)

# Start plotting
dotplot_df2 = dplyr::filter(dotplot_df, eRegulon_name %in% tf_list)
dotplot_df2[['eRegulon_name']] <- factor(dotplot_df2[['eRegulon_name']], levels = tf_list) 
base_plot <- ggplot(dotplot_df2, aes_string(y = 'celltype', x = 'eRegulon_name'))
base_plot <- base_plot +
  geom_tile(aes_string(fill = 'extended_gene_based_AUC')) +
  geom_point(aes_string(size = "extended_region_based_AUC"), color = "black") +
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  #scale_fill_carto_c(type = "diverging", palette = "TealRose", direction = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_blank())+
  coord_fixed()

ggsave(filename = "Reproducibility/Results/Plots/B/FigureS6E.pdf", plot = base_plot, width = 15, height = 7)