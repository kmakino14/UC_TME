#####################################################
# DOGMA analysis - MSC -
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
lineage = 'MSC'

DOGMA = file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")) %>% readRDS()

#################################################
# Fig.3B  Module violin plot & heatmap
#################################################
# Heatmap
DefaultAssay(DOGMA) = "RNA"

gene_list = list(
    MP1_genes = c("PDGFRA","LUM","PDPN","MMP2","TSHZ2"),
    MP2_genes = c("RGS5","ACTA2","TAGLN","PDGFA","MCAM"),
    MP3_genes = c("DES","MYH11","CASQ2","ADRA1A","DTNA"),
    MP4_genes = c("C3","C7","CFD","PI16","CXCL12"),
    MP5_genes = c("FAP","COL1A1","COL3A1","MMP14","TWIST1")
    )

genes.use = unlist(gene_list)

pseudobulk_tmp = AverageExpression(DOGMA, assay = "RNA", slot = "data", group.by = "celltype", 
                                   features = genes.use, return.seurat = TRUE)
pseudobulk_tmp$celltype = colnames(pseudobulk_tmp) %>% 
                          factor(., levels = c("proCAF","iCAF_CD321","iCAF_SLC14A1","matCAF","matCAP","contCAP","vSMC"))
pseudobulk_tmp = ScaleData(pseudobulk_tmp, scale.max = 10)

plots <- list()
for(tmp_num in 1:5){
p1 = DoHeatmap(pseudobulk_tmp, assay = "RNA", features = rownames(pseudobulk_tmp)[{5*(tmp_num)-4}:{5*(tmp_num)}], 
               group.by = "celltype", draw.lines = FALSE,
               slot = "scale.data", group.bar = TRUE, label = TRUE, 
               size = 3, hjust = 0, angle = 90, group.bar.height = 0.02, 
               disp.max = 2, disp.min = -2,
               group.colors=c("proCAF"       = "#1f77b4",
                              "iCAF_CD321"   = "#ff7f0e",
                              "iCAF_SLC14A1" = "#2ca02c",
                              "matCAF"       = "#d62728",
                              "matCAP5"      = "#9467bd",
                              "contCAP"      = "#8c564b",
                              "vSMC"         = "#e377c2")
               ) + 
     scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
                          midpoint = 0, guide = "colourbar", aesthetics = "fill",limits=c(-2.5,2.5)) + coord_fixed() + NoLegend() + theme(text = element_text(size = 20))

plots[[tmp_num]] = p1
}

#-------------------------------
# Violin plot

Hotspot_df = fread_n("Reproducibility/Results/Hotspot/MSC/UC_DOGMA_MSC_module_scores.txt") %>% t() %>% as.matrix()
DOGMA[["Hotspot"]]   = CreateAssayObject(data = Hotspot_df)

DefaultAssay(DOGMA) = "Hotspot"
DOGMA$celltype = factor(DOGMA$celltype, levels = c("proCAF","iCAF_CD321","iCAF_SLC14A1","matCAF","matCAP","contCAP","vSMC"))
cols = c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2")
p1  = VlnPlot(DOGMA, c("MP1"), group.by = "celltype", pt.size = 0, cols = cols) + NoLegend() 
p2  = VlnPlot(DOGMA, c("MP2"), group.by = "celltype", pt.size = 0, cols = cols) + NoLegend() 
p3  = VlnPlot(DOGMA, c("MP3"), group.by = "celltype", pt.size = 0, cols = cols) + NoLegend() 
p4  = VlnPlot(DOGMA, c("MP4"), group.by = "celltype", pt.size = 0, cols = cols) + NoLegend() 
p5  = VlnPlot(DOGMA, c("MP5"), group.by = "celltype", pt.size = 0, cols = cols) + NoLegend() 

combined_vlnplots <- plot_grid(p1, p2, p3, p4, p5, ncol = 5)

common_theme_vln <- theme(
  plot.margin = unit(c(0, 0, 0, 1), "cm"),
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
)

common_theme_hp <- theme(
  plot.margin = unit(c(0, 0, 0, 1), "cm"),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
)

# Apply the common theme to the violin plots and heatmaps
p1_vln <- p1 + common_theme_vln 
p2_vln <- p2 + common_theme_vln 
p3_vln <- p3 + common_theme_vln 
p4_vln <- p4 + common_theme_vln 
p5_vln <- p5 + common_theme_vln 

p1_heatmap <- plots[[1]] + common_theme_hp
p2_heatmap <- plots[[2]] + common_theme_hp
p3_heatmap <- plots[[3]] + common_theme_hp
p4_heatmap <- plots[[4]] + common_theme_hp
p5_heatmap <- plots[[5]] + common_theme_hp

final_combined_figure = plot_grid(p1_vln,p2_vln, p3_vln, p4_vln, p5_vln,
                                  p1_heatmap, p2_heatmap, p3_heatmap, p4_heatmap, p5_heatmap,
                                  nrow = 2, align = 'hv', rel_widths = c(1, 1, 1, 1, 1), rel_heights = c(0.4, 1))
# Save the final combined figure
ggsave("Reproducibility/Results/Plots/MSC/Figure3B.pdf", plot = final_combined_figure, width = 15, height = 10)


#################################################
# Fig.S4B  Module corrplot
#################################################
sig_df_1 = fread_n("Reproducibility/Results/VISIONR/UC_DOGMA_MSC_signature_score_hotspot.txt") %>% as.data.frame()
sig_df_2 = fread_n("Reproducibility/Results/VISIONR/UC_DOGMA_MSC_signature_score_literature.txt") %>% as.data.frame()
sig_df = cbind(sig_df_1, sig_df_2)

## corr
corr = cor(sig_df ,method="pearson")
order_tmp = hclust_col_order(corr, METHOD="ward.D2")
corr = corr[order_tmp, order_tmp]

palette_colors <- jdb_palette("brewer_jamaica")
color_scale <- colorRamp2(c(-1,-0.7,-0.4,-0.1,0.2,0.5,0.7,0.85,1), palette_colors[seq(length(palette_colors))])
#color_scale <- colorRamp2(seq(-1, 1, length.out = length(palette_colors)), palette_colors)

## Main Heatmap
p = Heatmap(corr, 
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

paste0("Reproducibility/Results/Plots/MSC/FigureS4B.pdf") %>% pdf(.,width = 5, height = 5)
 p  
dev.off()


#################################################
# Fig.3G  Milo beeswarm plot
#################################################
res = fread_n("Reproducibility/Results/Milo/output/Milo_MSC_design_Organ_STAGE2_contrasts_MI_vs_NMI_result.txt")

p = res %>%
    dplyr::filter(nhood_annotation_frac > 0.5) %>%
    mutate(signif=ifelse(SpatialFDR < 0.1, logFC, 0)) %>%
    group_by(nhood_annotation)%>%
    mutate(mean_lfc_val = ifelse(!is.na(signif), logFC, 0)) %>%
    mutate(mean_lfc = mean(mean_lfc_val)) %>%
    ungroup() %>%
    dplyr::arrange(mean_lfc) %>%
    mutate(nhood_annotation=factor(nhood_annotation, levels=c('vSMC',"contCAP",'matCAP','matCAF','iCAF_SLC14A1','iCAF_CD321','proCAF'))) %>%
    ggplot(aes(nhood_annotation, logFC, color=signif)) + 
        ggbeeswarm::geom_quasirandom(size=2) +
            coord_flip() +
            scale_color_gradient2(high='#8E063B', mid='grey', low='#023FA5', 
                                  name='DA logFC\n(10% SpatialFDR)',
                                  limits = c(-7, 7)) +
            theme_bw(base_size=18) +
            geom_hline(yintercept=0, linetype=2) +
            xlab('') + 
            ylim(-7, 7) 

paste0("Reproducibility/Results/Plots/MSC/Figure3G.pdf") %>% pdf(.,width = 7, height = 3.5)
 p  
dev.off()


#################################################
# Fig.3H  TF activity rank plot
#################################################
# 1) matCAF
tmp_result = paste0("Reproducibility/Results/LINGER/Primary/cell_population_exp_TF_activity_zscore_matCAF.txt") %>%
             fread_n() %>%
             dplyr::arrange(., desc(TF)) %>%
             rownames_to_column("TF_name") %>% 
             rownames_to_column("rank")
tmp_result$rank = as.numeric(tmp_result$rank)

label = c('TWIST1','RUNX1','PRRX1','SMAD2')
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
labs(title = "matCAF", y = "activity score", x = "Rank Sorted Annotations")

pdf('Reproducibility/Results/Plots/MSC/Figure3H_matCAF.pdf', width = 5, height = 3)
 print(p1)
dev.off()

# 2) matCAP
tmp_result = paste0("Reproducibility/Results/LINGER/Primary/cell_population_exp_TF_activity_zscore_matCAP.txt") %>%
             fread_n() %>%
             dplyr::arrange(., desc(TF)) %>%
             rownames_to_column("TF_name") %>% 
             rownames_to_column("rank")
tmp_result$rank = as.numeric(tmp_result$rank)

label = c('TWIST1','RUNX1','PRRX1','SMAD2')
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
labs(title = "matCAP", y = "activity score", x = "Rank Sorted Annotations")

pdf('Reproducibility/Results/Plots/MSC/Figure3H_matCAP.pdf', width = 5, height = 3)
 print(p2)
dev.off()


#########################################################
## Fig.3I/S4F  Gene track with co-accessible regions
#########################################################
source("Reproducibility/Scripts/Source/Seurat_source.R")

celltype_list=c("proCAF","iCAF_CD321","iCAF_SLC14A1","matCAF","matCAP","contCAP","vSMC")

region_list <- list(
  COL3A1 = c('chr2-188900000-189200000'),
  COL5A2 = c('chr2-188900000-189200000'),
  CXCL14 = c('chr5-135542000-135580000'),
  CD34 = c("chr1-207880000-207925000")
)

motif_list = list(
    COL3A1 = c('SMAD2','RUNX1','PRRX1'),
    COL5A2 = c('SMAD2','RUNX1','PRRX1'),
    CXCL14 = c('PITX1','FOS'),
    CD34 = c('KLF4')
)

cols <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2")


Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='MSC',group_list=celltype_list,gene='COL3A1',
               region_list=region_list, motif_list=motif_list,cutoff=0.3,cols=cols,
               path='Reproducibility/Results/Plots/MSC/FIgure3I_COL3A1.pdf')

Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='MSC',group_list=celltype_list,gene='COL5A2',
               region_list=region_list, motif_list=motif_list,cutoff=0.3,cols=cols,
               path='Reproducibility/Results/Plots/MSC/FIgure3I_COL5A2.pdf')

Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='MSC',group_list=celltype_list,gene='CXCL14',
               region_list=region_list, motif_list=motif_list,cutoff=0.2,cols=cols,
               path='Reproducibility/Results/Plots/MSC/FIgureS4F_CXCL14.pdf')

Plot_genetrack(Obj=DOGMA,group.by='celltype',lineage='MSC',group_list=celltype_list,gene='CD34',
               region_list=region_list, motif_list=motif_list,cutoff=0.1,cols=cols,
               path='Reproducibility/Results/Plots/MSC/FIgureS4F_CD34.pdf')


#########################################################
## Fig.S4D  proCAF Kaplan-Meier
#########################################################

TCGA_surv = fread_n("Reproducibility/Data/External_cohort/TCGA_BLCA_metadata.txt")
result = fread_n("Reproducibility/Results/BayesPrism/bayesprism_TCGA_prop_df.txt")

result$pct_proCAF = result$proCAF/(result$proCAF+result$iCAF)*100
result$patient = paste0(take_factor(rownames(result),1,'-'),'-',take_factor(rownames(result),2,'-'),'-',
                        take_factor(rownames(result),3,'-'),'-',take_factor(rownames(result),4,'-'))
result <- result[grepl("01A$", result$patient), ]
result$patient <- substr(result$patient, 1, nchar(result$patient) - 1)
result <- result[!duplicated(result$patient), ]
rownames(result) = result$patient

result = result[rownames(TCGA_surv),]

df <- cbind(result, TCGA_surv)

res.cut <- tryCatch({
    surv_cutpoint(df, time = "OS_time", event = "OS_event", minprop = 0.1, variables = 'pct_proCAF')
}, error = function(e) {
    NULL  # Return NULL if an error occurs
})

if (is.null(res.cut)) {
    next  # Skip the rest of the loop and move to the next iteration
}

res.cat <- surv_categorize(res.cut)

fit <- survfit(as.formula(paste("Surv(OS_time, OS_event) ~", tmp_module)), data = res.cat)
cox_fit <- coxph(as.formula(paste("Surv(OS_time, OS_event) ~", tmp_module)), data = res.cat)
summary_cox <- summary(cox_fit)
    
p_value <- summary_cox$sctest["pvalue"]
hazard_ratio <- 1 / summary_cox$coefficients[1, "exp(coef)"]
formatted_p_value <- format.pval(p_value, digits = 2, eps = .Machine$double.eps)

# Create the survival plot with formatted p-value
p1 <- ggsurvplot(
    fit, 
    data = res.cat, 
    risk.table = TRUE, 
    conf.int = TRUE,
    pval = paste("p =", formatted_p_value),   # Adds the formatted p-value to the plot
    pval.method = FALSE,  # Optionally remove the test method label if not needed
    palette = c("#ef7c21", "#2372a9"),  # Use custom colors for the survival curves
    pval.coord = c(max(fit$time) * 0.5, 1), # Customize the position of the p-value
    censor.shape = 124,
    censor.size = 2  
)
p1$plot <- p1$plot + 
  annotate("text", x = max(fit$time) * 0.75, y = 0.8, 
           label = paste("HR =", round(hazard_ratio, 3)), 
           size = 5)

paste0('Reproducibility/Results/Plots/MSC/FIgureS4D.pdf') %>% pdf_3(., w = 6, h = 6)
  print(p1)
dev.off()