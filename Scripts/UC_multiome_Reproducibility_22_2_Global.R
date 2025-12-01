#####################################################
# DOGMA analysis - Global -
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

#################################################
# QC violinplot by coarse celltype
#################################################
metadata = fread_n("Reproducibility/Data/UC_DOGMA_metadata.txt")
metadata$coarse_celltype2 = fct_collapse(metadata$coarse_celltype,
    Epithelial = c('UC','NEC'),
    MSC = c('CAF','Pericyte')
    )
metadata$coarse_celltype2 = factor(metadata$coarse_celltype2,
    levels = rev(c("CD4_Tconv",'Treg','CD8_T','NK_ILC','B','Mono_Mac',"DC","Mast","Epithelial",'MSC','Endothelial'))
    )

metadata["log10_UMI_RNA"] = log10(metadata["total_counts_RNA"])
metadata["log10_Genes_RNA"] = log10(metadata["n_genes_RNA"])
metadata["log10_total_fragment_counts"] = log10(metadata["total_counts_ATAC"])

p1 = ggplot(metadata, aes(x = coarse_celltype2, y = log10_UMI_RNA, fill = coarse_celltype2))+
  geom_violin(trim = FALSE)+
  geom_boxplot(width = .1, color = "white", outlier.color = NA)+
  scale_fill_manual(values=rev(c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#1B9E77","#D95F02","#E7298A"))) +
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

p2 = ggplot(metadata, aes(x = coarse_celltype2, y = log10_Genes_RNA, fill = coarse_celltype2))+
  geom_violin(trim = FALSE)+
  geom_boxplot(width = .1, color = "white", outlier.color = NA)+
  scale_fill_manual(values=rev(c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#1B9E77","#D95F02","#E7298A"))) +
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

p3 = ggplot(metadata, aes(x = coarse_celltype2, y = pct_counts_mt, fill = coarse_celltype2))+
  geom_violin(trim = FALSE)+
  geom_boxplot(width = .1, color = "white", outlier.color = NA)+
  scale_fill_manual(values=rev(c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#1B9E77","#D95F02","#E7298A"))) +
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

p4 = ggplot(metadata, aes(x = coarse_celltype2, y = log10_total_fragment_counts, fill = coarse_celltype2))+
  geom_violin(trim = FALSE)+
  geom_boxplot(width = .1, color = "white", outlier.color = NA)+
  scale_fill_manual(values=rev(c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#1B9E77","#D95F02","#E7298A"))) +
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()

p <- plot_grid(p1,p2,p3,p4,ncol=4)

paste0("Reproducibility/Results/Plots/Global/FigureS1F.pdf") %>% pdf(., w = 14, h = 4)
 plot(p)
dev.off()


#################################################
# proportion boxplot by STAGE
#################################################
metadata <- metadata %>% mutate(coarse_celltype3 = if_else(celltype == "Normal", "Normal", coarse_celltype2))
metadata$coarse_celltype3 = fct_collapse(metadata$coarse_celltype3, Malignant = c('Epithelial'))

metadata_IM = dplyr::filter(metadata, coarse_celltype3 %in% c("CD4_Tconv","Treg","CD8_T","NK_ILC","B","Mono_Mac","DC","Mast"))
metadata_IM$coarse_celltype3 = metadata_IM$coarse_celltype3[,drop=TRUE]
metadata_IM$coarse_celltype3 = factor(metadata_IM$coarse_celltype3, levels = rev(c("CD4_Tconv","Treg","CD8_T","NK_ILC","B","Mono_Mac","DC","Mast")))

metadata_NI = dplyr::filter(metadata, coarse_celltype3 %in% c("Normal","Malignant","MSC","Endothelial"))
metadata_NI$coarse_celltype3 = metadata_NI$coarse_celltype3[,drop=TRUE] 
metadata_NI$coarse_celltype3 = factor(metadata_NI$coarse_celltype3, levels = rev(c("Normal","Malignant","MSC","Endothelial")))

################################
metadata_IM = dplyr::select(metadata_IM, c('coarse_celltype3','sample','STAGE'))
ct_levels <- levels(metadata_IM$coarse_celltype3)
sample_stage <- metadata_IM %>% dplyr::distinct(sample, STAGE)

cell_frac <- metadata_IM %>%
  dplyr::count(sample, coarse_celltype3, name = "n") %>%
  left_join(sample_stage, by = "sample") %>% 
  group_by(sample, STAGE) %>%
  complete(coarse_celltype3 = ct_levels, fill = list(n = 0)) %>%
  mutate(prop = 100*n / sum(n)) %>% 
  ungroup() %>%
  dplyr::select(sample, coarse_celltype3, STAGE, prop)

cell_frac$coarse_celltype3 = factor(cell_frac$coarse_celltype3, levels = c("CD4_Tconv","Treg","CD8_T","NK_ILC","B","Mono_Mac","DC","Mast"))
cell_frac$STAGE            = factor(cell_frac$STAGE,            levels = c("Early_L","Early_H","Advanced","post_BCG"))

stage_comparisons <- combn(levels(cell_frac$STAGE), 2, simplify = FALSE)
palette <- c("#00AFBB", "#E7B800", "#FC4E07",'#CC61B0')

# Create boxplot grouped by STAGE_v2 for each variable, with pairwise comparisons
p = ggplot(cell_frac, aes(x = STAGE, y = prop, fill = STAGE)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6, shape = 16) +
  scale_fill_manual(values = palette) +
  facet_wrap(~ coarse_celltype3, scales = "free_y", ncol = 8) +
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

################################
BCG_pt = c("BC_027","BC_032","BC_039","BC_040","BC_043","BC_044","BC_047","BC_048")
metadata_NI = dplyr::filter(metadata_NI, !sample %in% BCG_pt) 
metadata_NI$sample = metadata_NI$sample[,drop=TRUE]

metadata_NI = dplyr::select(metadata_NI, c('coarse_celltype3','sample','STAGE'))
ct_levels <- levels(metadata_NI$coarse_celltype3)
sample_stage <- metadata_NI %>% dplyr::distinct(sample, STAGE)

cell_frac <- metadata_NI %>%
  dplyr::count(sample, coarse_celltype3, name = "n") %>%
  left_join(sample_stage, by = "sample") %>% 
  group_by(sample, STAGE) %>%
  complete(coarse_celltype3 = ct_levels, fill = list(n = 0)) %>%
  mutate(prop = 100*n / sum(n)) %>% 
  ungroup() %>%
  dplyr::select(sample, coarse_celltype3, STAGE, prop)

cell_frac$coarse_celltype3 = factor(cell_frac$coarse_celltype3, levels = c("Normal","Malignant","MSC","Endothelial"))
cell_frac$STAGE            = factor(cell_frac$STAGE,            levels = c("Early_L","Early_H","Advanced"))

stage_comparisons <- combn(levels(cell_frac$STAGE), 2, simplify = FALSE)
palette2 <- c("#00AFBB", "#E7B800", "#FC4E07")

# Create boxplot grouped by STAGE_v2 for each variable, with pairwise comparisons
q = ggplot(cell_frac, aes(x = STAGE, y = prop, fill = STAGE)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6, shape = 16) +
  scale_fill_manual(values = palette2) +
  facet_wrap(~ coarse_celltype3, scales = "free_y", ncol = 4) +
  labs(x = "Status", y = "Propotion") +
  theme_classic() +
  theme(legend.position = "none") +
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

combined_plot <- p + q +
  plot_layout(widths = c(1, 0.502))

paste0("Reproducibility/Results/Plots/Global/FigureS1G.pdf") %>% pdf(., w = 24.8, h = 3)
 combined_plot
dev.off()
