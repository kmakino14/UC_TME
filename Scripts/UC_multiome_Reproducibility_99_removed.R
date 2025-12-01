#################################################
# Fig.S2F  Module usage by STAGE
#################################################
hs_df = fread_n("Reproducibility/Results/Hotspot/Malignant/UC_DOGMA_Malignant_module_scores.txt")
metadata = fread_n("Reproducibility/Data/UC_DOGMA_metadata.txt") %>% dplyr::filter(., celltype %in% c("LUM",'NRP','SQM','MES','NEC'))
small_number_samples = c("BC_011","BC_012","BC_014","BC_030","BC_033","BC_042","BC_045")

hs_use = c("MP1","MP2","MP3","MP4","MP5","MP7","MP8","MP9","MP10","MP11")
hs_df_w_meta = cbind(hs_df, metadata) %>% 
      dplyr::select(.,c(hs_use, 'sample', 'STAGE')) %>% 
      dplyr::filter(., !sample %in% small_number_samples)

hs_df_avg <- hs_df_w_meta %>%
  group_by(sample) %>%
  summarise(across(-STAGE, mean, na.rm = TRUE)) %>% as.data.frame()

# Find unique, consistent relationships between 'sample' and 'STAGE'
consistent_relationships <- hs_df_w_meta %>%
  group_by(sample) %>%
  dplyr::filter(n_distinct(STAGE) == 1) %>%
  distinct(sample, STAGE) %>% as.data.frame()

combined_df <- hs_df_avg %>%
  left_join(consistent_relationships, by = "sample")

# Convert data to long format for easy plotting
combined_long <- combined_df %>%
  pivot_longer(cols = c(MP1,MP2,MP3,MP4,MP5,MP7,MP8,MP9,MP10,MP11), names_to = "Variable", values_to = "Value")

# Define the factor levels for Variable to set the facet order
combined_long$Variable <- factor(combined_long$Variable, levels = hs_use)

# Define the factor levels for STAGE and color palette
combined_long$STAGE <- factor(combined_long$STAGE, levels = c('Early_L', 'Early_H', 'Advanced'))
stage_comparisons <- combn(levels(combined_long$STAGE), 2, simplify = FALSE)
palette <- c("#00AFBB", "#E7B800", "#FC4E07")

# Create boxplot grouped by STAGE for each variable, with pairwise comparisons
p = ggplot(combined_long, aes(x = STAGE, y = Value, fill = STAGE)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6, shape = 16) +
  scale_fill_manual(values = palette) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 11) +
  labs(x = "Stage", y = "Value") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) +
  geom_pwc(
    aes(group = STAGE_v2),
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "fdr",
    bracket.nudge.y = 0.1,
    hide.ns = TRUE,
    tip.length = 0
  )

paste0("Reproducibility/Results/Plots/Malignant/FigureS2F.pdf") %>% pdf_3(., w=22, h=2.5)
 plot(p)
dev.off()

#--------------------------------------------------------------------------------------------

DefaultAssay(DOGMA) = 'RNA'

# Load necessary libraries
df <- data.frame(
Group = DOGMA@meta.data$celltype,
Value = DOGMA[["RNA"]]@data['F3',]
)

# Create the plot
df$Group = factor(df$Group, levels = c('LUM','NRP','SQM','MES','NEC'))
colors <- c("LUM" = "#533f7c", "NRP" = "#37888e", "SQM" = "#4c8d49", "MES" = "#ce8d24", "NEC" = "#aa4743")

# Create the plot with violin, boxplot, and pairwise comparison
p <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
     geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with fill by group
     geom_boxplot(width = 0.14, fill = "white", outlier.shape = NA) +  # Boxplot overlay
     scale_fill_manual(values = colors) +  # Apply custom colors
     theme_classic() +
     labs(title = 'F3', y = paste0("Normalized CD142 expression")) +
     theme(legend.position = "none")

p = p + stat_compare_means(
    comparisons = list(c("LUM", "NRP"),c("SQM", "NRP"),c("MES", "NRP"),c("NEC", "NRP")),
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "fdr",
    bracket.nudge.y = 0.1,
    hide.ns = TRUE,
    tip.length = 0
    )

paste0('Reproducibility/Results/Plots/Malignant/FigureS3F_RNA.pdf') %>% pdf(.,h=4,w=5)
 plot(p)
dev.off()

#--------------------------------------------------------------------------------------------

adata.obs["dummy_all_nan"] = np.nan

# Plot — all points will be colored "#F2F2F2"
sc.pl.embedding(
    adata,
    basis="umap",
    color=["dummy_all_nan"],
    color_map=cmap,
    show=False
)

plt.savefig(f"{fig_dir}FigureS6D_UMAP_pseudotime_bg.pdf", bbox_inches='tight')
plt.close()

↓ bugになったら以下を使用

# Use a constant value instead of NaN
adata.obs["dummy_all_nan"] = 0.0

from matplotlib.colors import ListedColormap
import scanpy as sc

# Colormap with a single light gray color
cmap_grey = ListedColormap(["#F2F2F2"])

sc.pl.embedding(
    adata,
    basis="umap",
    color="dummy_all_nan",   # note: string, not list
    color_map=cmap_grey,
    vmin=0,
    vmax=0,
    show=True,
)

#--------------------------------------------------------------------------------------------

adata.obs["dummy_all_nan"] = np.nan

# Plot — all points will be colored "#F2F2F2"
sc.pl.embedding(
    adata,
    basis="umap",
    color=["dummy_all_nan"],
    color_map=cmap,
    show=False
)

plt.savefig(f"{fig_dir}Figure4D_UMAP_pseudotime_bg.pdf", bbox_inches='tight')
plt.close()

#--------------------------------------------------------------------------------------------

df <- read.csv("Reproducibility/Results/CellRank2/CD4_T/fate_probabilities_dpt_Treg.csv", row.names = 1)
prob_score = df[, c("Treg_effector_2", "Treg_effector_3")] 
prob_score2 = prob_score-0.5 

# Create Seurat object with dummy counts
seurat_obj <- CreateSeuratObject(counts = prob_score2 %>% t(), meta.data = df)

# Attach UMAP coordinates
seurat_obj[["UMAP"]] <- CreateDimReducObject(
  embeddings = df[, c("UMAP_1", "UMAP_2")] %>% as.matrix(),
  key = "umap_",
  assay = "RNA"
)

p = SCpubr::do_FeaturePlot(
  seurat_obj,
  features = "Treg-effector-2",
  enforce_symmetry = TRUE,
  diverging.palette = "RdYlGn",
  diverging.direction = -1,
  plot_cell_borders = TRUE,
  pt.size = 4,
  border.size = 1,
  raster = TRUE,
  raster.dpi = 1024)

paste0("Reproducibility/Results/Plots/CD4_T/Figure4D_fate_prob_UMAP.pdf") %>% pdf(w=3, h=5)
 plot(p)
dev.off()


#--------------------------------------------------------------------------------------------
#########################################################
## Fig.S7X  Gene track with co-accessible regions
#########################################################
source("Reproducibility/Scripts/Source/Seurat_source.R")

metadata = fread_n('Reproducibility/Results/Milo/output/Milo_Treg_effector_metadata.txt')
DOGMA_Treg = subset(DOGMA, cells = rownames(metadata))
DOGMA_Treg$nhood_groups = metadata[colnames(DOGMA_Treg),]$nhood_groups

Idents(DOGMA_Treg) = 'nhood_groups'
DOGMA_part = subset(DOGMA_Treg, idents = c('in_nhoods_NMI_enr','in_nhoods_MI_enr'))
DOGMA_part$nhood_groups = factor(DOGMA_part$nhood_groups, levels = c('in_nhoods_NMI_enr','in_nhoods_MI_enr'))

group_list=c('in_nhoods_NMI_enr','in_nhoods_MI_enr')
region_list <- list(
  TBX21 = c('chr17-47697000-47750000')
)

motif_list = list(
  TBX21 = c('TBX21','BATF')
)

cols <- c('#D6ED17FF', '#3E282BFF')

Plot_genetrack_nhood(Obj=DOGMA_part,lineage='CD4_T',group_list=group_list,gene='TBX21',
                     region_list=region_list, motif_list=motif_list,cutoff=0.15,cols=cols,
                     path='Reproducibility/Results/Plots/CD4_T/FIgureS7X.pdf')


#--------------------------------------------------------------------------------------------
# P06
colors <- c('P06_stress_Hi' = "#b12e23",
           'P06_clone_1' = '#c2b32e',
           'P06_clone_2' = '#bae1e3',
           'P06_clone_3' = '#eeada8',
           'P06_clone_4' = "#ea9223")

plots <- list()
for(tmp_TF in c("KLF5_extended_+/+_(113g)","FOSL2_extended_+/+_(122g)")){
   tmp_TF_core = paste0(take_factor(tmp_TF,1,"_"), '_', take_factor(tmp_TF,2,"_"))
   plot_df <- eRegulon_df %>%
     dplyr::mutate(celltype = eRegulon_meta$clone3) %>%
     dplyr::filter(., celltype %in% c('P06_stress_Hi','P06_clone_1','P06_clone_2','P06_clone_3','P06_clone_4'))
   plot_df$celltype = factor(plot_df$celltype, levels = c('P06_stress_Hi','P06_clone_1','P06_clone_2','P06_clone_3','P06_clone_4'))
   
   plots[[tmp_TF]] = ggplot(plot_df, aes(x = celltype, y = .data[[tmp_TF]])) +
     geom_boxplot(fill = "white", outlier.shape = NA) +  # suppress outliers
     geom_quasirandom(
       aes(color = celltype),
       method = 'pseudorandom',
       size = 3,
       alpha = 0.75,
       shape = 16
     ) +
     labs(
       title = paste0(tmp_TF),
       x = "Cell Type",
       y = "eRegulon score"
     ) +  scale_color_manual(values = colors) +
     theme_classic() +
     theme(legend.position = "none")
     # stat_compare_means(comparisons = comparisons, label = "p.signif", p.adjust.method = 'BH')
}

paste0('Reproducibility/Results/Plots/Slide-tags/FigureS13E_beeswarm_plot.pdf') %>% pdf_3(w=5, h=3.5)
 plots
dev.off()

#--------------------------------------------------------------------------------------------

palette_clone2 = ['#d3d3d3', '#f9d6a5', '#91cb98', '#d3d3d3', '#d3d3d3', '#d3d3d3', '#d3d3d3', '#d3d3d3', 
                  '#d3d3d3', '#d3d3d3', '#d3d3d3', "#d3d3d3", '#d3d3d3', '#d3d3d3', '#d3d3d3', "#d3d3d3",
                  '#d3d3d3', '#d3d3d3', '#d3d3d3', '#d3d3d3', '#005084']

labels = adata_sp.obs['lineage_w_clone'].cat.categories.tolist()
color_dict = dict(zip(labels, palette_clone2))

pucks = adata_sp_L.copy()
tdatas = {sample: pucks[df.index] for sample, df in pucks.obs.groupby('specimen') }
axs = state_plot_grid_L(padding=1.5)
fig = tc.pl.scatter(tdatas,
                    keys="lineage_w_clone",
                    colors=color_dict,
                    position_key=['x','y'],
                    joint=True, point_size=4, 
                    ax=axs, 
                    noticks=False,
                    legend=True,
                    rasterized=False,
                    render=False
                    )
plt.savefig(f"Reproducibility/Results/Plots/Slide-tags/Figure7E.pdf", bbox_inches='tight')
plt.close()