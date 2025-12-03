#####################################################
# TREKKER Seurat object construction
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
suppressMessages(library(SCpubr))

suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsci))
suppressMessages(library(colorspace))
suppressMessages(library(BuenColors))

suppressMessages(library(survival))
suppressMessages(library(survminer))

suppressMessages(library(dplyr))
set.seed(1234)

#******************************************
data_dir="Reproducibility/Data"

matrix_files = paste0(data_dir, '/matrix_files/TREKKER/GEX/')
meta = fread_n(paste0(data_dir, "/TREKKER/UC_TREKKER_metadata.txt"))

RNA_counts <- Read10X(data.dir = matrix_files)
scRNA = CreateSeuratObject(counts = RNA_counts , assay = "RNA")
scRNA = NormalizeData(scRNA, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
scRNA = AddMetaData(scRNA, metadata = meta)
scRNA = ScaleData(scRNA, assay = "RNA")

latent_df = fread_n(paste0(data_dir, "/embeddings/UC_TREKKER_Global_MultiVI_latent_df.txt"))
scRNA[['scvi_latent']] <- CreateDimReducObject(embeddings = as.matrix(latent_df), key = "scvi_", assay = "RNA")
umap_df = fread_n(paste0(data_dir, "/embeddings/UC_TREKKER_Global_MultiVI_UMAP_df.txt"))
scRNA[['umap']] <- CreateDimReducObject(embeddings = as.matrix(umap_df), key = "umap_", assay = "RNA")

saveRDS(scRNA, paste0(data_dir,"/Seurat/UC_TREKKER_seurat_obj_Global_RNA.rds"))

#################################################
# VISIONR
#################################################

scRNA = paste0(data_dir,"/Seurat/UC_TREKKER_seurat_obj_Global_RNA.rds") %>% readRDS()

Idents(scRNA) = 'celltype'
scRNA = subset(scRNA, ident =c('LUM','SQM','NRP'))

TREKKER_signature_score_literature <- function(scRNA){
	DefaultAssay(scRNA) = 'RNA'
	gene_df = fread(paste0("Reproducibility/Data/Signature_gene_list/TableS3_module_genes_literature_Malignant.txt")) %>% as.data.frame()
	
	# Convert data.frame / data.table to list of columns
	sig_list <- lapply(gene_df, function(x) {
	  x <- as.character(x)
	  x <- x[x != ""]        # remove empty strings
	  x <- na.omit(x)        # remove NA
	  unique(x)
	})
	
	signatures=c()
	for(tmp_num in 1:length(sig_list)){
	    tmp_genes = sig_list[[tmp_num]]
	    tmp_genes2 = intersect(tmp_genes, rownames(scRNA))
	    tmp_module = c(rep(1, length(tmp_genes2)))
	    names(tmp_module) = tmp_genes2
	    sig <- createGeneSignature(name = names(sig_list)[[tmp_num]], sigData = tmp_module)
	    signatures = c(signatures,sig)
	}
	  
	vision.obj <- Vision(scRNA,
	                     signatures = signatures,
	                     dimRed = "scvi_latent",
	                     projection_methods = c("UMAP"),
	                     meta = scRNA@meta.data)
	
	options(mc.cores = 4)
	vis <- analyze(vision.obj)
	
	sig_df = getSignatureScores(vis) %>% as.data.frame()
	write.table_n_2(sig_df, 'CB', paste0('Reproducibility/Results/VISIONR/UC_TREKKER_Malignant_signature_score_literature.txt'))
}

TREKKER_signature_score_hotspot <- function(scRNA){
	DefaultAssay(scRNA) = 'RNA'

    gene_df = fread(paste0("Reproducibility/Data/Signature_gene_list/TableS3_module_genes_Malignant_hotspot.txt")) %>% as.data.frame()
    gene_df <- gene_df[gene_df$Annotation != "" & !is.na(gene_df$Annotation), ]
  
    sig_list <- gene_df %>%
    group_by(Annotation) %>%
    summarise(genes = list(Gene)) %>%
    deframe()
  
    names(sig_list) <- gsub(" ", "_", names(sig_list))
  
    signatures=c()
    for(tmp_num in 1:length(sig_list)){
        tmp_genes = sig_list[[tmp_num]]
        tmp_genes2 = intersect(tmp_genes, rownames(scRNA))
        tmp_module = c(rep(1, length(tmp_genes2)))
        names(tmp_module) = tmp_genes2
        sig <- createGeneSignature(name = names(sig_list)[[tmp_num]], sigData = tmp_module)
        signatures = c(signatures,sig)
    }
    
    vision.obj <- Vision(scRNA,
                         signatures = signatures,
                         dimRed = "scvi_latent",
                         projection_methods = c("UMAP"),
                         meta = scRNA@meta.data)
    
    options(mc.cores = 4)
    vis <- analyze(vision.obj)
    
    sig_df = getSignatureScores(vis) %>% as.data.frame()
	write.table_n_2(sig_df, 'CB', paste0('Reproducibility/Results/VISIONR/UC_TREKKER_Malignant_signature_score_hotspot.txt'))
}

TREKKER_signature_score_literature(scRNA)
TREKKER_signature_score_hotspot(scRNA)

#################################################
# Fig.7H/7I/S13F　Violin plots 
#################################################
sig_df1 = fread_n('Reproducibility/Results/VISIONR/UC_TREKKER_Malignant_signature_score_literature.txt')
sig_df2 = fread_n('Reproducibility/Results/VISIONR/UC_TREKKER_Malignant_signature_score_hotspot.txt')

sig_df = cbind(sig_df1, sig_df2)
sig_df = sig_df[colnames(scRNA),]
scRNA@meta.data = cbind(scRNA@meta.data, sig_df)

#-------------------------------------------
plots <- list()
modules = c('Hypoxia','Stress_Barkley','Partial_epithelial-mesenchymal_transition','Interferon_signaling')

for(tmp_module in modules){
	df <- data.frame(
	  Group = scRNA@meta.data$clone2,
	  Value = scRNA@meta.data[tmp_module]
	)
	colnames(df) = c('Group','Value')
	df = dplyr::filter(df, Group %in% c('P02_clone_1','P02_clone_1_hypoxic','P02_clone_2'))
	df$Group = factor(df$Group, levels = c('P02_clone_1_hypoxic','P02_clone_1','P02_clone_2'))
	colors <- c('P02_clone_1_hypoxic' = "#c64f1b", 'P02_clone_1' = '#f4d2a3', 'P02_clone_2' = '#8dc494')
	
	# Create the plot with violin, boxplot, and pairwise comparison
	plots[[tmp_module]] <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
	  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with fill by group
	  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +  # Boxplot overlay
	  scale_fill_manual(values = colors) +  # Apply custom colors
	  theme_classic() +
	  labs(title = tmp_module, y = paste0("Signature score")) +
	  theme(legend.position = "none")
	
	my_comparisons <- list(c('P02_clone_1_hypoxic', 'P02_clone_1'),c('P02_clone_1_hypoxic', 'P02_clone_2'),c('P02_clone_1', 'P02_clone_2'))
	
	plots[[tmp_module]] = plots[[tmp_module]] + stat_compare_means(
	    comparisons = my_comparisons,
	    method = "wilcox.test",
	    label = "p.signif",
	    p.adjust.method = "BH",
	    bracket.nudge.y = 0.1,
	    hide.ns = TRUE,
	    tip.length = 0
	  )
}

paste0("Reproducibility/Results/Plots/Slide-tags/Figure7H.pdf") %>% pdf(.,h=4,w=5)
 plots
dev.off()

#-------------------------------------------
plots <- list()
modules = c('Stress_Barkley')

for(tmp_module in modules){
	df <- data.frame(
	  Group = scRNA@meta.data$clone2,
	  Value = scRNA@meta.data[tmp_module]
	)
	colnames(df) = c('Group','Value')
	df = dplyr::filter(df, Group %in% c('P06_clone_1','P06_clone_1_stress','P06_clone_2','P06_clone_2_stress','P06_clone_3','P06_clone_4'))
	df$Group = factor(df$Group, levels = c('P06_clone_1_stress','P06_clone_1','P06_clone_2_stress','P06_clone_2','P06_clone_3','P06_clone_4'))
	colors <- c('P06_clone_1_stress' = "#f71804",
	            'P06_clone_1' = '#d2bf1b',
	            'P06_clone_2_stress' ="#b12e23",
	            'P06_clone_2' = '#BFA4D6',
	            'P06_clone_3' = '#F4B1AB',
	            'P06_clone_4' = "#94B0CB")
	
	# Create the plot with violin, boxplot, and pairwise comparison
	plots[[tmp_module]] <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
	  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with fill by group
	  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +  # Boxplot overlay
	  scale_fill_manual(values = colors) +  # Apply custom colors
	  theme_classic() +
	  labs(title = tmp_module, y = paste0("Signature score")) +
	  theme(legend.position = "none")
}

paste0("Reproducibility/Results/Plots/Slide-tags/FigureS13F_signature_scores.pdf") %>% pdf(.,h=4,w=5)
 plots
dev.off()

#-------------------------------------------
Idents(scRNA) = 'sample'
scRNA_P06 = subset(scRNA, ident = c('P06'))
Idents(scRNA_P06) = 'celltype'
scRNA_P06 = subset(scRNA_P06, ident = c('LUM'))
scRNA_P06$clone2 = factor(scRNA_P06$clone2, levels = c('P06_clone_1_stress','P06_clone_1','P06_clone_2_stress','P06_clone_2','P06_clone_3','P06_clone_4'))
colors <- c('P06_clone_1_stress' = "#f71804",
            'P06_clone_1' = '#d2bf1b',
            'P06_clone_2_stress' ="#b12e23",
            'P06_clone_2' = '#BFA4D6',
            'P06_clone_3' = '#F4B1AB',
            'P06_clone_4' = "#94B0CB")

p1 <- SCpubr::do_ViolinPlot(sample = scRNA_P06, 
                           assay = 'RNA',
                           slot = 'data',
                           group.by = 'clone2',
                           line_width = 0,
                           features = 'DUSP1',
                           plot_boxplot = FALSE,
                           order = FALSE,
                           colors.use = colors)
p2 <- SCpubr::do_ViolinPlot(sample = scRNA_P06, 
                           assay = 'RNA',
                           slot = 'data',
                           group.by = 'clone2',
                           line_width = 0,
                           features = 'FOS',
                           plot_boxplot = FALSE,
                           order = FALSE,
                           colors.use = colors)

paste0("Reproducibility/Results/Plots/Slide-tags/FigureS13F_RNA.pdf") %>% pdf(.,h=8,w=5)
 plot(p1/p2)
dev.off()

#-------------------------------------------
TF_activity = fread_n("Reproducibility/Results/LINGER/TREKKER/output/cell_population_TF_activity.txt") %>% t() %>% as.data.frame()
TF_activity = TF_activity[colnames(scRNA),]

plots <- list()
modules = c('NFE2L2', 'HIF1A')

for(tmp_module in modules){
  # Load necessary libraries
  df <- data.frame(
    Group = scRNA@meta.data$clone2,
    Value = TF_activity[tmp_module]
  )
  colnames(df) = c('Group','Value')
  df = dplyr::filter(df, Group %in% c('P02_clone_1','P02_clone_1_hypoxic','P02_clone_2'))
  
  # Create the plot
  df$Group = factor(df$Group, levels = c('P02_clone_1_hypoxic','P02_clone_1','P02_clone_2'))
  colors <- c(
        	'P02_clone_1_hypoxic' = "#c64f1b",
        	'P02_clone_1' = '#f4d2a3',
        	'P02_clone_2' = '#8dc494')
  
  # Create the plot with violin, boxplot, and pairwise comparison
  plots[[tmp_module]] <- ggplot(df, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with fill by group
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +  # Boxplot overlay
    scale_fill_manual(values = colors) +  # Apply custom colors
    theme_classic() +
    labs(title = tmp_module, y = paste0("TF activity")) +
    theme(legend.position = "none")
  
  # Add pairwise comparisons for specific groups using geom_pwc()
  comparisons <- list(c('P02_clone_1_hypoxic', 'P02_clone_1'),c('P02_clone_1_hypoxic', 'P02_clone_2'),c('P02_clone_1', 'P02_clone_2'))
  
  plots[[tmp_module]] = plots[[tmp_module]] + stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif",
      p.adjust.method = "BH",
      bracket.nudge.y = 0.1,
      hide.ns = TRUE,
      tip.length = 0
    )
}

paste0('Reproducibility/Results/Plots/Slide-tags/Figure7I.pdf') %>% pdf(.,h=4,w=5)
 plots
dev.off()


#################################################
# Fig.7J　Beeswarm plots of eRegulons
#################################################
eRegulon_df = fread_n("Reproducibility/Results/scenicplus/TREKKER/outs/eRegulon_df.txt") # metacell x TF 
eRegulon_meta = fread_n("Reproducibility/Results/scenicplus/TREKKER/outs/eRegulon_metadata.txt")
colnames(eRegulon_meta) = c('clone2', 'clone3','ATAC_clone2', 'ATAC_clone3')
cols_with_plus <- grep("\\+/\\+", colnames(eRegulon_df), value = TRUE)

# P02
plots <- list()
for(tmp_TF in c("NFE2L2_extended_+/+_(376g)")){
    tmp_TF_core = paste0(take_factor(tmp_TF,1,"_"), '_', take_factor(tmp_TF,2,"_"))
    plot_df <- eRegulon_df %>%
      dplyr::mutate(celltype = eRegulon_meta$clone3) %>%
      dplyr::filter(., celltype %in% c('P02_clone_1_hypoxic','P02_clone_1','P02_clone_2'))
    plot_df$celltype = factor(plot_df$celltype, levels = c('P02_clone_1_hypoxic','P02_clone_1','P02_clone_2'))
    
    comparisons <- list(c('P02_clone_1_hypoxic','P02_clone_1'), c('P02_clone_1','P02_clone_2'), 
    					c('P02_clone_1_hypoxic','P02_clone_2'))
    
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
        x = "Clone",
        y = "eRegulon score"
      ) +  scale_color_manual(values = c(
      	'P02_clone_1_hypoxic' = "#c64f1b",
      	'P02_clone_1' = '#f4d2a3',
      	'P02_clone_2' = '#8dc494')) +
      theme_classic() +
      theme(legend.position = "none") +
      stat_compare_means(comparisons = comparisons, label = "p.signif", p.adjust.method = 'BH')
}

paste0('Reproducibility/Results/Plots/Slide-tags/Figure7J.pdf') %>% pdf_3(w=5, h=3.5)
 plots
dev.off()

#################################################
# Fig.7C  Scatter plot of TF expression/activity
#################################################
TF_exp_df = fread_n("Reproducibility/Results/LINGER/TREKKER/output/cell_population_exp_TF_activity_zscore_NRP.txt") %>% 
						rownames_to_column('Name')
						
highlight_pos <- c("SNAI2", "FOXP2", "POU6F2", "NFIA", "RFX3", 'LHX9', 'TWIST2', 'ARNTL', 'CREM', 'TCF7L1',
                   "MEF2C", 'YY1', 'STAT3')

highlight_neg <- c('GATA3', 'PPARG', 'FOXA1')

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

paste0('Reproducibility/Results/Plots/Slide-tags/Figure7C.pdf') %>% pdf_3(w=3, h=3)
  p
dev.off()


#################################################
# Fig.7G/S13G　Co-occurence heatmap
#################################################
# P02
df = fread_n("Reproducibility/Results/TREKKER/co_occurrence_matrix_P02_zscore_distance_50.txt")

col_list = c("#6889D0",'#A3B4E5',"#D3DBF4", "#F6F6F6", 
           "#F7D3D3", "#E6A4A4", "#CB6F70", "#9D3D3D")

path_pheatmap = paste0("Reproducibility/Results/Plots/Slide-tags/Figure7G.pdf")
col = colorRampPalette(col_list)
col_fun <- colorRamp2(seq(-1, 1, length.out = 50), col(50))


pheatmap::pheatmap(df,
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   cellwidth = 10, cellheight = 10,
                   color = col(50),
                   breaks = seq(-20,20,length.out=50),
                   filename=path_pheatmap,
                   border_color = "black" ,na_col = "gray90",
                   main = "co_occurrence_P02_zscore_distance_50")

###################################################
# P06
df = fread_n("Reproducibility/Results/TREKKER/co_occurrence_matrix_P06_zscore_distance_50.txt")

col_list = c("#6889D0",'#A3B4E5',"#D3DBF4", "#F6F6F6", 
           "#F7D3D3", "#E6A4A4", "#CB6F70", "#9D3D3D")

path_pheatmap = paste0("Reproducibility/Results/Plots/Slide-tags/FigureS13G.pdf")
col = colorRampPalette(col_list)
col_fun <- colorRamp2(seq(-1, 1, length.out = 50), col(50))

pheatmap::pheatmap(df,
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   cellwidth = 10, cellheight = 10,
                   color = col(50),
                   breaks = seq(-5,5,length.out=50),
                   filename=path_pheatmap,
                   border_color = "black" ,na_col = "gray90",
                   main = "co_occurrence_P06_zscore_distance_50")
