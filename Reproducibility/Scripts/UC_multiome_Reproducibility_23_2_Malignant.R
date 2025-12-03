#####################################################
# DOGMA analysis - Malignant -
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
lineage = 'Malignant'

DOGMA = file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")) %>% readRDS()

#################################################
# Fig.2C  Module heatmap
#################################################
hs_df = fread_n("Reproducibility/Results/Hotspot/Malignant/UC_DOGMA_Malignant_module_scores.txt")
metadata = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_metadata.txt") %>% dplyr::filter(., celltype %in% c("LUM",'NRP','SQM','MES','NEC'))

hs_df_part = hs_df[,c("MP1","MP2","MP3","MP4","MP5")]
colnames(hs_df_part) = c("LUM",'NRP',"SQM","MES","NEC")

# Normalize each module score
norm = function(data){scale(data, center = min(data), scale = (max(data) - min(data)))}
hs_df_part_norm = apply(hs_df_part,2,norm)
rownames(hs_df_part_norm) = rownames(hs_df_part)

norm_part = hs_df_part_norm %>% t()
#zscore_part[zscore_part>2]  = 2
#zscore_part[-2>zscore_part] = -2

metadata$celltype = factor(metadata$celltype, levels=c("LUM",'NRP',"SQM","MES","NEC"))
group_list        = c("LUM",'NRP',"SQM","MES","NEC")

# reoderding in each celltype
CB_order = c()
for(i in 1:length(group_list)){
    tmp_data = norm_part[,rownames(metadata)[metadata$celltype==group_list[i]]] %>% t() %>% as.data.frame()
    if(i==1){
        tmp_data_reorder = dplyr::arrange(tmp_data, desc(LUM)) %>% rownames()
    }else if(i==2){
        tmp_data_reorder = dplyr::arrange(tmp_data, desc(NRP)) %>% rownames()
    }else if(i==3){
        tmp_data_reorder = dplyr::arrange(tmp_data, desc(SQM)) %>% rownames()
    }else if(i==4){
        tmp_data_reorder = dplyr::arrange(tmp_data, desc(MES)) %>% rownames()
    }else{
        tmp_data_reorder = dplyr::arrange(tmp_data, desc(NEC)) %>% rownames()
    }
    CB_order = c(CB_order, tmp_data_reorder)
}

metadata2 = metadata[CB_order,]
data = norm_part[c('LUM','NRP','SQM','MES','NEC'),CB_order]

# Create Heatmap
metadata2$celltype <- factor(metadata2$celltype, levels = c('LUM', 'NRP', 'SQM', 'MES', 'NEC'))
celltype_colors <- c("LUM" = "#533f7c", "NRP" = "#37888e", "SQM" = "#4c8d49", "MES" = "#ce8d24", "NEC" = '#aa4743')

p2 <- HeatmapAnnotation(df = metadata2[, "celltype", drop = FALSE],
                        col = list(celltype = celltype_colors),
                        simple_anno_size = unit(0.4, "cm"),
                        annotation_name_gp = gpar(fontsize = 8))

col_fun = circlize::colorRamp2(c(0, 0.5, 1), plasma(99)[c(1, 60, 99)])

p1 = Heatmap(data, 
             col = col_fun,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             row_names_gp = gpar(fontsize = 3),
             column_names_gp = gpar(fontsize = 0),
             top_annotation = p2,
             heatmap_legend_param = list(title = "Hotspot", at = c(0, 0.5, 1)),
             use_raster = TRUE
             )


paste0("Reproducibility/Results/Plots/Malignant/Figure2C_heatmap.pdf") %>% pdf(.,h=2, w=5)
 p1  
dev.off()


#################################################
# Fig.S2C  Module corrplot
#################################################
sig_df_1 = fread_n("Reproducibility/Results/VISIONR/UC_DOGMA_Malignant_signature_score_hotspot.txt") %>% as.data.frame()
sig_df_2 = fread_n("Reproducibility/Results/VISIONR/UC_DOGMA_Malignant_signature_score_literature.txt") %>% as.data.frame()
colnames(sig_df_1) = paste0(colnames(sig_df_1), '_de_novo')
sig_df = cbind(sig_df_1, sig_df_2)

#============================
# Lineage program 
#============================
sig_df_lineage = dplyr::select(sig_df, c('Luminal_de_novo','Neural-like_progenitor_de_novo',
    'Neuroendocrine_de_novo','Mesenchymal_de_novo','Squamous_de_novo',
    'NEC','LUM', 'NRP', 'MES','Mesenchymal','SQM','BSL'))

## corr
corr = cor(sig_df_lineage ,method="pearson")
order_tmp = hclust_col_order(corr, METHOD="ward.D2")
corr = corr[order_tmp, order_tmp]

palette_colors <- jdb_palette("brewer_jamaica")
color_scale <- colorRamp2(seq(-1, 1, length.out = length(palette_colors)), palette_colors)

## Main Heatmap
p1 = Heatmap(corr, 
             col = color_scale, 
             rect_gp = gpar(type = "none"),
             cluster_rows = FALSE,                    ## clustering dendrogram
             cluster_columns = FALSE,
             cell_fun = function(i, j, x, y, w, h, fill) {
                 if(i <= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                }
             },
             row_names_gp = gpar(fontsize = 6),
             column_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Pearson r", at = c(-1, 0, 1)),
             heatmap_width = unit(8, "cm"), heatmap_height = unit(8, "cm"),
             use_raster = TRUE
             )

paste0("Reproducibility/Results/Plots/Malignant/FigureS2C_lineage.pdf") %>% pdf(.,width = 5, height = 5)
 p1  
dev.off()

#============================
# State program 
#============================
sig_df_state = dplyr::select(sig_df, c('Cycling/G2M_de_novo','Interferon_signaling_de_novo',
                'Stress_de_novo','Partial_epithelial-mesenchymal_transition_de_novo',
                'Cycle','CYG',"EMT-I",'Stress_Gavish','Stress_Barkley',
                "Interferon_MHC-II-I",'Interferon','MP7'))

## corr
corr = cor(sig_df_state ,method="pearson")
order_tmp = hclust_col_order(corr, METHOD="ward.D2")
corr = corr[order_tmp, order_tmp]

palette_colors <- jdb_palette("brewer_celsius")
color_scale <- colorRamp2(seq(-1, 1, length.out = length(palette_colors)), palette_colors)

## Main Heatmap
p2 = Heatmap(corr, 
             col = color_scale, 
             rect_gp = gpar(type = "none"),
             cluster_rows = FALSE,                    ## clustering dendrogram
             cluster_columns = FALSE,
             cell_fun = function(i, j, x, y, w, h, fill) {
                 if(i <= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                }
             },
             row_names_gp = gpar(fontsize = 6),
             column_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Pearson r", at = c(-1, 0, 1)),
             heatmap_width = unit(8, "cm"), heatmap_height = unit(8, "cm"),
             use_raster = TRUE
             )

paste0("Reproducibility/Results/Plots/Malignant/FigureS2C_state.pdf") %>% pdf(.,width = 5, height = 5)
 p2  
dev.off()


#################################################
# Fig.S2F  Module usage by Origin
#################################################
hs_df = fread_n("Reproducibility/Results/Hotspot/Malignant/UC_DOGMA_Malignant_module_scores.txt")
metadata = fread_n("Reproducibility/Data/DOGMA/UC_DOGMA_metadata.txt") %>% dplyr::filter(., celltype %in% c("LUM",'NRP','SQM','MES','NEC'))
small_number_samples = c("BC_011","BC_012","BC_014","BC_030","BC_033","BC_042","BC_045")
hs_use = c("MP6","MP7","MP8","MP9","MP10","MP11")

hs_df_w_meta = cbind(hs_df, metadata) %>% 
      dplyr::select(.,c(hs_use,'sample', 'Organ')) %>% 
      dplyr::filter(., !sample %in% small_number_samples)

hs_df_avg <- hs_df_w_meta %>%
  group_by(sample) %>%
  summarise(across(-Organ, mean, na.rm = TRUE)) %>% as.data.frame()

# Find unique, consistent relationships between 'sample' and 'Organ'
consistent_relationships <- hs_df_w_meta %>%
  group_by(sample) %>%
  dplyr::filter(n_distinct(Organ) == 1) %>%
  distinct(sample, Organ) %>% as.data.frame()

combined_df <- hs_df_avg %>%
  left_join(consistent_relationships, by = "sample")

# Convert data to long format for easy plotting
combined_long <- combined_df %>%
  pivot_longer(cols = c(MP6,MP7,MP8,MP9,MP10,MP11), names_to = "Variable", values_to = "Value")

# Define the factor levels for Variable to set the facet order
combined_long$Variable <- factor(combined_long$Variable, levels = hs_use)

combined_long$Organ <- factor(combined_long$Organ, levels = c('BC', 'UTUC'))
stage_comparisons <- combn(levels(combined_long$Organ), 2, simplify = FALSE)
palette <- c("#A8E6A3", "#C5A3C6")

# Create boxplot
p = ggplot(combined_long, aes(x = Organ, y = Value, fill = Organ)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6, shape = 16) +
  scale_fill_manual(values = palette) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 11) +
  labs(x = "Organ", y = "Value") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) +
  geom_pwc(
    aes(group = Organ),
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "fdr",
    bracket.nudge.y = 0.1,
    hide.ns = TRUE,
    tip.length = 0
  )

paste0("Reproducibility/Results/Plots/Malignant/FigureS2F.pdf") %>% pdf_3(., w=10, h=2.5)
 plot(p)
dev.off()


#################################################
# Fig.S2G  Module usage by celltype
#################################################
## import module score calculated by Hotspot
hs_df = fread_n("Reproducibility/Results/Hotspot/Malignant/UC_DOGMA_Malignant_module_scores.txt") %>%
        dplyr::select(., c("MP7","MP8","MP9","MP10","MP11"))

hs_df_shifted <- as.data.frame(lapply(hs_df, function(column) {
  min_val <- min(column, na.rm = TRUE)
  shifted_column <- column - min_val + 1
  return(shifted_column)
}))
rownames(hs_df_shifted) = rownames(hs_df)

DOGMA[["Hotspot_shifted"]]   = CreateAssayObject(data = hs_df_shifted %>% as.matrix() %>% t()) ## "_" -> "-" に変換される

cols = rev(c('#533f7c','#37888e','#4c8d49','#ce8d24','#aa4743'))
DOGMA$celltype = factor(DOGMA$celltype, levels = rev(c("LUM","NRP","SQM","MES","NEC")))
p = RidgePlot(object = DOGMA, assay = 'Hotspot_shifted', cols = cols, 
              features = c("MP7","MP8","MP9","MP10","MP11"), 
              group.by = 'celltype', log = TRUE)

paste0("Reproducibility/Results/Plots/Malignant/FigureS2G.pdf") %>% pdf(., w = 12, h = 6)
 plot(p)
dev.off()


#################################################
# Fig.2D  Marker peaks heatmap
#################################################
marker_peaks <- fromJSON('Reproducibility/Results/SnapATAC2/Malignant/snapatac2_Malignant_marker_peaks_filtered_p_0.05.json')
peak_mat_z = fread_n("Reproducibility/Results/SnapATAC2/Malignant/snapatac2_Malignant_filtered_marker_peak_zscored.txt")   # peaks*zscored_log2(1+RPKM)

new_order = c('LUM', 'NRP', 'SQM', 'MES', 'NEC')
peak_mat_z = peak_mat_z[,new_order]

set.seed(123)
n_sample = 50000
sampled_df <- peak_mat_z[sample(nrow(peak_mat_z), n_sample), ]

# Initialize an empty list to store the reordered dataframes
reordered_dfs <- list()

# Loop through each element in marker_peaks
for (column_name in new_order) {
  target_rows <- marker_peaks[[column_name]]
  target_rows = intersect(target_rows,rownames(sampled_df))
  sub_df <- sampled_df[target_rows, , drop = FALSE]
  sub_df <- sub_df[order(sub_df[[column_name]], decreasing = TRUE), ]
  reordered_dfs[[column_name]] <- sub_df
}

# Concatenate all the reordered dataframes
final_df <- do.call(rbind, reordered_dfs)

# Define a custom RdBu palette with lighter blue and unchanged red
light_rd_bu_palette <- colorRampPalette(c("#88cce7", "#e0f3f8", "white", "#fdae61", "#d73027"))(99)

# Generate the heatmap with the modified palette
p  = Heatmap(final_df[rev(rownames(final_df)),], 
             col = light_rd_bu_palette,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             row_names_gp = gpar(fontsize = 0),
             column_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "log2(1+RPKM)", at = c(-4, 0, 4)),
             use_raster = FALSE
             )

# Save as PDF
paste0("Reproducibility/Results/Plots/Malignant/Figure2D.pdf") %>% pdf(., width = 4, height = 4)
 draw(p)
dev.off()


#################################################
# Fig.S3B  Pancancer neuroendcrine signature
#################################################
hs_df = fread_n("Reproducibility/Results/VISIONR/UC_DOGMA_Malignant_signature_score_literature.txt") %>%
        as.data.frame() %>% dplyr::select(., c("Subtype-N","Subtype-P"))

DOGMA[["NEC"]] = CreateAssayObject(data = hs_df %>% t() %>% as.matrix()) 

subtypes = c("Subtype-N" ,"Subtype-P")

group_levels <- c("LUM","NRP","SQM","MES","NEC")
colors <- c("LUM"="#533f7c","NRP"="#37888e","SQM"="#4c8d49","MES"="#ce8d24","NEC"="#aa4743")

# Build one long data frame for all signatures
df_all <- purrr::map_dfr(subtypes, function(sig) {
  tibble(
    Group = DOGMA@meta.data$celltype,
    Value = as.numeric(DOGMA[["NEC"]]@data[sig, , drop = TRUE]),
    Signature = sig
  )
})

df_all <- df_all %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Signature = factor(Signature, levels = subtypes)
  )

# One plot with facets (one panel per signature)
p <- ggplot(df_all, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.14, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  facet_wrap(~ Signature, scales = "free_y", ncol = 3) +
  theme_classic() +
  labs(y = "Normalized expression", x = NULL) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold")
  )

# Save once
ggsave("Reproducibility/Results/Plots/Malignant/FigureS3B.pdf", p, width = 6, height = 2.5)


#################################################
# Fig.S3A  TF activity rank plot
#################################################
# 1) MES
tmp_result = paste0("Reproducibility/Results/LINGER/Primary/cell_population_exp_TF_activity_zscore_MES.txt") %>%
             fread_n() %>%
             dplyr::arrange(., desc(TF)) %>%
             rownames_to_column("TF_name") %>% 
             rownames_to_column("rank")
tmp_result$rank = as.numeric(tmp_result$rank)

label = c('PAX3','ETV1','SOX10','ASCL1','PHOX2B','NEUROD1','POU3F3','ZIC2','NFIB','TCF4')
tmp_result$label <- ifelse(tmp_result$TF_name %in% label, tmp_result$TF_name, NA)
labeled_points <- tmp_result[!is.na(tmp_result$label) & tmp_result$label != "", ]

custom_palette <- jdb_palette("brewer_spectra")
color_range <- range(tmp_result$TF)
max_abs_sum <- max(abs(color_range))  # Find the maximum absolute value

p1 <- ggplot(tmp_result, aes(x = rank, y = TF, label = label, color = TF)) +
scale_x_continuous(breaks=c(1,200, 400, 600),limits=c(-150,650)) + 
scale_y_continuous(breaks=seq(-2.5,4,by=0.5)) + 
geom_point(size = 0.5, shape = 16, alpha = 0.5) + 
geom_point(data = labeled_points, aes(x = rank, y = TF), size = 1, shape = 21, color="black") + 
scale_color_gradientn(colors = custom_palette, 
                      limits = c(-max_abs_sum, max_abs_sum), 
                      values = scales::rescale(c(-max_abs_sum, 0, max_abs_sum))) +
geom_text_repel(color = "black",
                size = 1,
                segment.color = "black",
                segment.size = 0.25,
                min.segment.length = 0,
                na.rm = TRUE,
                show.legend = TRUE) +
theme_classic() + 
labs(title = "MES", y = "activity score", x = "Rank Sorted Annotations")

pdf('Reproducibility/Results/Plots/Malignant/FigureS3A_MES.pdf', width = 5, height = 3)
 print(p1)
dev.off()

# 2) NEC
tmp_result = paste0("Reproducibility/Results/LINGER/Primary/cell_population_exp_TF_activity_zscore_NEC.txt") %>%
             fread_n() %>%
             dplyr::arrange(., desc(TF)) %>%
             rownames_to_column("TF_name") %>% 
             rownames_to_column("rank")
tmp_result$rank = as.numeric(tmp_result$rank)

label = c('POU2F3','PROX1','SOX4','E2F1','TFAP2B','SOX2')
tmp_result$label <- ifelse(tmp_result$TF_name %in% label, tmp_result$TF_name, NA)
labeled_points <- tmp_result[!is.na(tmp_result$label) & tmp_result$label != "", ]

custom_palette <- jdb_palette("brewer_spectra")
color_range <- range(tmp_result$TF)
max_abs_sum <- max(abs(color_range))  # Find the maximum absolute value

p2 <- ggplot(tmp_result, aes(x = rank, y = TF, label = label, color = TF)) +
scale_x_continuous(breaks=c(1,200, 400, 600),limits=c(-150,650)) + 
scale_y_continuous(breaks=seq(-2.5,4,by=0.5)) + 
geom_point(size = 0.5, shape = 16, alpha = 0.5) + 
geom_point(data = labeled_points, aes(x = rank, y = TF), size = 1, shape = 21, color="black") +
scale_color_gradientn(colors = custom_palette, 
                      limits = c(-max_abs_sum, max_abs_sum), 
                      values = scales::rescale(c(-max_abs_sum, 0, max_abs_sum))) +
geom_text_repel(color = "black",
                size = 1,
                segment.color = "black",
                segment.size = 0.25,
                min.segment.length = 0,
                na.rm = TRUE,
                show.legend = TRUE) +
theme_classic() + 
labs(title = "NEC", y = "activity score", x = "Rank Sorted Annotations")

pdf('Reproducibility/Results/Plots/Malignant/FigureS3A_NEC.pdf', width = 5, height = 3)
 print(p2)
dev.off()


#################################################
# Fig.2F/Fig.S3E  Multicox regression
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

Get_signature_score_z = function(study,module_df){
  if(study == 'TCGA'){
    counts_tmp = fread_n("Reproducibility/Data/External_cohort/TCGA_BLCA_RNAseq_data.txt")
    metadata = fread_n("Reproducibility/Data/External_cohort/TCGA_BLCA_metadata.txt")
  } else{
    counts_tmp = fread_n("Reproducibility/Data/External_cohort/IMvigor210_RNAseq_data.txt")
    metadata = fread_n("Reproducibility/Data/External_cohort/IMvigor210_metadata.txt")
  }

  ## ssGSEA
  counts = t(counts_tmp)
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
  df_ssGSEA_w_meta = cbind(df_ssGSEA_z,metadata)
  return(df_ssGSEA_w_meta)
}

TCGA_ssGSEA_z    = Get_signature_score_z(study = 'TCGA', module_df = module_df)
IMvigor_ssGSEA_z = Get_signature_score_z(study = 'IMvigor', module_df = module_df)

#================
# TCGA
#================
df_for_cox = dplyr::select(TCGA_ssGSEA_z, -PFS_event, -PFS_time) %>%
             transmute(OS_time,
               OS_event,
               Age,
               Sex = factor(Sex, labels = c("Male", "Female")),
               STAGE = factor(STAGE, labels = c("Stage II","Stage III","Stage IV")),
               LUM,NRP,SQM,MES,NEC,MTN,CYG,IFN,SEC,STR,pEMT
               )
plot = forest_model(coxph(Surv(OS_time, OS_event) ~ ., df_for_cox))

paste0('Reproducibility/Results/Plots/Malignant/Figure2F.pdf') %>% pdf_3(., w=10,h=10)
  plot
dev.off()

#================
# IMvigor210
#================
df_for_cox = dplyr::select(IMvigor_ssGSEA_z, -Received_platinum, -Lund2) %>%
             transmute(os,
                       censOS,
                       Sex,
                       Baseline_ECOG_Score,
                       Tissue,
                       LUM,NRP,SQM,MES,NEC,MTN,CYG,IFN,SEC,STR,pEMT)
plot = forest_model(coxph(Surv(os, censOS) ~ ., df_for_cox))

paste0('Reproducibility/Results/Plots/Malignant/FigureS3E.pdf') %>% pdf_3(., w=10,h=10)
  plot
dev.off()

#################################################
# Fig.2G  Volcano plot (NRP vs LUM)
#################################################

data = fread_n('Reproducibility/Results/Differential/UC_DOGMA_DEG_MAST_Malignant_NRP_vs_LUM.txt')

label_up = c("WNT5B",'BRINP3','LRFN5','FOXP2','F3','ROBO1', 'CSPG4')
label_down = c("GATA3","FOXA1","PPARG","UPK1B","KRT7","CDH1")

paste0('Reproducibility/Results/Plots/Malignant/Figure2G.pdf') %>% pdf(., w=3,h=3)
  Volcano_DEG_label(data = data, celltype = "NRP", label_up = label_up, label_down = label_down, log2FC = 0.5)
dev.off()


#################################################
# Fig.S2F  CD142 violon plot
#################################################
DefaultAssay(DOGMA) = 'ADT'

# Load necessary libraries
df <- data.frame(
Group = DOGMA@meta.data$celltype,
Value = DOGMA[["ADT"]]@data['surface-A0822-CD142',]
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
     labs(title = 'CD142', y = paste0("Normalized CD142 expression")) +
     theme(legend.position = "none")

p = p + stat_compare_means(
    comparisons = list(c("LUM", "NRP"),c("SQM", "NRP"),c("MES", "NRP"),c("NEC", "NRP")),
    method = "wilcox.test",
    label = "p.signif",
    p.adjust.method = "BH",
    bracket.nudge.y = 0.1,
    hide.ns = TRUE,
    tip.length = 0
    )

paste0('Reproducibility/Results/Plots/Malignant/FigureS3F.pdf') %>% pdf(.,h=4,w=5)
 plot(p)
dev.off()


#################################################
# Fig.S3I  TF activity/expression scatter
#################################################
TF_exp_df = fread_n('Reproducibility/Results/LINGER/Primary/cell_population_exp_TF_activity_zscore_NRP.txt') %>%
            rownames_to_column('Name')

highlight_pos <- c("SNAI2", "FOXP2", "POU6F2","HNF4G", "NFIA", 'TWIST2', 
                   "MEF2C", "STAT3", "RFX1", 'YY1')
highlight_neg <- c('GATA3', 'POU2F3', 'NPAS2')

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
     theme(legend.position = "none",
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.line = element_line(color = "black"),
           panel.border = element_rect(color = "black", fill = NA))

paste0("Reproducibility/Results/Plots/Malignant/FigureS3I.pdf") %>% pdf(., w=3, h=3)
 plot(p)
dev.off()


#########################################################
## Fig.2J/S3H  Kaplan-Meier by TF activity
#########################################################
library(survival)
library(survminer)
library(purrr)

## Load & align
TF_activity <- t(fread_n("Reproducibility/Results/LINGER/Primary/TCGA_BLCA_TF_activity_general.txt"))
TCGA_surv   <- fread_n("Reproducibility/Data/External_cohort/TCGA_BLCA_RNAseq_FPKM_for_LINGER_metadata.txt")

TF_activity <- TF_activity[rownames(TCGA_surv), , drop = FALSE]
TCGA_data   <- cbind(TCGA_surv, TF_activity)

## Cutpoint & categorize (for ALL TFs at once)
res.cut <- surv_cutpoint(
  TCGA_data,
  time = "OS_MONTHS",
  event = "OS_STATUS",
  variables = colnames(TF_activity)
)
res.cat_all <- surv_categorize(res.cut)                       # only categorized TFs
res.cat <- cbind(TCGA_surv[, c("OS_MONTHS","OS_STATUS")],     # put back survival columns
                 res.cat_all[, colnames(TF_activity), drop=FALSE])

plot_km <- function(var, out_pdf) {
  df  <- res.cat[, c("OS_MONTHS","OS_STATUS", var)]
  colnames(df)[3] <- "grp"

  fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ grp, data = df)
  cox <- coxph   (Surv(OS_MONTHS, OS_STATUS) ~ grp, data = df)
  s   <- summary(cox)

  pval <- s$coefficients[1, "Pr(>|z|)"]
  hr   <- 1/s$coefficients[1, "exp(coef)"]

  p <- ggsurvplot(
    fit, data = df,
    risk.table = TRUE, conf.int = TRUE,
    palette = c("#ef7c21", "#2372a9"),
    pval = paste("p =", format.pval(pval, digits = 2)),
    pval.coord = c(max(fit$time) * 0.5, 1),
    censor.shape = 124, censor.size = 2
  )
  p$plot <- p$plot + annotate("text",
    x = max(fit$time) * 0.75, y = 0.8,
    label = paste("HR =", round(hr, 3)), size = 5
  )

  pdf_3(out_pdf, h = 5, w = 4.5); print(p); dev.off()
}

out_map <- c(
  FOXP2  = "Reproducibility/Results/Plots/Malignant/Figure2J.pdf",
  STAT3 = "Reproducibility/Results/Plots/Malignant/FigureS3J_STAT3.pdf",
  SNAI2  = "Reproducibility/Results/Plots/Malignant/FigureSJH_SNAI2.pdf"
)

tfs <- c("FOXP2", "STAT3", "SNAI2")
walk(tfs, ~ plot_km(.x, out_map[[.x]]))


#########################################################
## Fig.2I/S3C-D  Gene track with co-accessible regions
#########################################################
source("Reproducibility/Scripts/Source/Seurat_source.R")

celltype_list=c('LUM','NRP','SQM','MES','NEC')

region_list <- list(
  ROBO1 = c("chr3-78990000-79116000"),
  RGS13 = c("chr1-192620000-192665000"),
  PIK3R5 = c('chr17-8860000-8980000'),
  NEUROD1 = c("chr2-181650000-181690000"),
  PHOX2B = c("chr4-41737000-41800000")
)

motif_list = list(
    ROBO1 = c('FOXP2','POU6F2'),
    RGS13 = c('POU2F3'),
    PIK3R5 = c('PROX1','POU2F3'),
    NEUROD1 = c('NEUROD1','ASCL1','NFIB'),
    PHOX2B = c('NEUROD1','PHOX2B')
)

cols <- c('#564182', '#378D94', '#4C924B', '#D69121', '#B14845')


Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='Malignant',group_list=celltype_list,gene='NEUROD1',
               region_list=region_list, motif_list=motif_list,cutoff=0.1,cols=cols,
               path='Reproducibility/Results/Plots/Malignant/FIgureS3C_NEUROD1.pdf')

Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='Malignant',group_list=celltype_list,gene='PHOX2B',
               region_list=region_list, motif_list=motif_list,cutoff=0.1,cols=cols,
               path='Reproducibility/Results/Plots/Malignant/FIgureS3C_PHOX2B.pdf')

Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='Malignant',group_list=celltype_list,gene='PIK3R5',
               region_list=region_list, motif_list=motif_list,cutoff=0.1,cols=cols,
               path='Reproducibility/Results/Plots/Malignant/FIgureS3D_PIK3R5.pdf')

Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='Malignant',group_list=celltype_list,gene='RGS13',
               region_list=region_list, motif_list=motif_list,cutoff=0.35,cols=cols,
               path='Reproducibility/Results/Plots/Malignant/FIgureS3D_RGS13.pdf')

Plot_genetrack2(Obj=DOGMA,lineage='Malignant',celltype_list=celltype_list,gene='ROBO1',
               region_list=region_list, motif_list=motif_list,cutoff=0.3,cols=cols,
               path='Reproducibility/Results/Plots/Malignant/FIgure2I.pdf')
